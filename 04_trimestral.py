"""
Monsaraz em Mudança — passo 04
Serie TRIMESTRAL da area de agua, com mascaramento de nuvens.

Novidades face ao passo 03:
  - 4 trimestres por ano (2017-2025) => captura o ciclo sazonal.
  - Mascara nuvens/sombras pela banda SCL do Sentinel-2.
  - Combina as imagens do trimestre numa mediana (um dia nublado e
    "tapado" por outro dia limpo).
  - Reporta, por trimestre, a % da freguesia SEM dados validos (buraco
    deixado pelas nuvens) -> indicador de confianca.

Grava tudo em area_trimestral.csv.
"""

import ee
import json
import csv

ee.Initialize(project="projetobeta-513108")

# Limite da freguesia (CAOP)
with open("monsaraz.geojson", encoding="utf-8") as fh:
    geom = json.load(fh)["features"][0]["geometry"]
aoi = ee.Geometry(geom, None, False)

PROJ = "EPSG:3763"                       # PT-TM06, metrico
area_aoi_ha = aoi.area(1, ee.Projection(PROJ)).divide(10000).getInfo()
print("Area da freguesia: %.1f ha\n" % area_aoi_ha)

TRIMESTRES = {1: ("01-01", "03-31"), 2: ("04-01", "06-30"),
              3: ("07-01", "09-30"), 4: ("10-01", "12-31")}


def mascara_nuvens(img):
    """Remove nuvens, sombras de nuvem, cirrus e no-data pela banda SCL."""
    scl = img.select("SCL")
    # 0 no-data, 1 saturado, 3 sombra, 8/9 nuvem, 10 cirrus
    mau = (scl.eq(0).Or(scl.eq(1)).Or(scl.eq(3))
           .Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)))
    # ficar so com as bandas que usamos -> coleccao homogenea (B3, B11)
    return img.updateMask(mau.Not()).select(["B3", "B11"])


def agua_no_trimestre(ano, q):
    ini, fim = TRIMESTRES[q]
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
           .filterBounds(aoi)
           .filterDate(f"{ano}-{ini}", f"{ano}-{fim}")
           .map(mascara_nuvens))

    n = col.size()
    comp = col.median().clip(aoi)              # composicao mediana do trimestre

    mndwi = comp.normalizedDifference(["B3", "B11"]).rename("MNDWI")
    agua = mndwi.gt(0)

    pix = ee.Image.pixelArea()
    area_agua = (agua.selfMask().multiply(pix)
                 .reduceRegion(ee.Reducer.sum(), aoi, 10, PROJ, maxPixels=1e9)
                 .get("MNDWI"))
    # mascara de dados validos (1 onde a mediana tem dados)
    valido = comp.select("B3").mask()
    area_valida = (valido.multiply(pix)
                   .reduceRegion(ee.Reducer.sum(), aoi, 10, PROJ, maxPixels=1e9)
                   .get("B3"))

    info = ee.Dictionary({
        "n": n,
        "area_ha": ee.Number(area_agua).divide(10000),
        "valida_ha": ee.Number(area_valida).divide(10000),
    }).getInfo()
    gap = 100.0 * (1 - info["valida_ha"] / area_aoi_ha)
    return info["n"], info["area_ha"], gap


linhas = []
print("ano | trim | imgs | buraco | agua (ha)")
print("----+------+------+--------+----------")
for ano in range(2017, 2026):
    for q in (1, 2, 3, 4):
        try:
            n, area_ha, gap = agua_no_trimestre(ano, q)
            flag = "  <-- pouco fiavel" if gap > 20 else ""
            print("%4d|  T%d  | %4d | %5.1f%% | %8.1f%s" % (ano, q, n, gap, area_ha, flag))
            linhas.append([ano, q, n, round(gap, 1), round(area_ha, 1)])
        except Exception as e:
            print("%4d|  T%d  | erro: %s" % (ano, q, e))

with open("area_trimestral.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["ano", "trimestre", "n_imagens", "buraco_pct", "area_ha"])
    w.writerows(linhas)

print("\nGuardado: area_trimestral.csv (%d trimestres)" % len(linhas))
