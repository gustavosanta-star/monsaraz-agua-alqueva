"""
Monsaraz em Mudança — passo 07
INCERTEZA: sensibilidade da área ao limiar do MNDWI.

Para cada ano (T3), recalcula a área com o limiar em -0.05, 0 e +0.05.
A diferença entre os extremos dá uma banda de incerteza honesta para o indicador.
Baixar o limiar inclui mais pixéis (área maior); subir, menos.
"""

import ee
import json
import csv

ee.Initialize(project="projetobeta-513108")

with open("monsaraz.geojson", encoding="utf-8") as fh:
    geom = json.load(fh)["features"][0]["geometry"]
aoi = ee.Geometry(geom, None, False)

PROJ = "EPSG:3763"
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()
LIMIARES = [-0.05, 0.0, 0.05]


def mascara(img):
    scl = img.select("SCL")
    mau = (scl.eq(0).Or(scl.eq(1)).Or(scl.eq(3)).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)))
    return img.updateMask(mau.Not()).select(["B3", "B11"])


def mndwi_verao(ano):
    comp = (S2.filterBounds(aoi).filterDate(f"{ano}-07-01", f"{ano}-09-30")
            .map(mascara).median().clip(aoi))
    return comp.normalizedDifference(["B3", "B11"])


def area_ha(mask):
    m2 = (mask.selfMask().multiply(pix)
          .reduceRegion(ee.Reducer.sum(), aoi, 10, PROJ, maxPixels=1e9).values().get(0))
    return ee.Number(ee.Algorithms.If(m2, m2, 0)).divide(10000)


linhas = []
print("ano | -0.05 (ha) | 0 (ha) | +0.05 (ha) | banda ± (ha)")
print("----+-----------+--------+------------+-------------")
for ano in range(2017, 2026):
    mndwi = mndwi_verao(ano)
    d = ee.Dictionary({
        "baixo": area_ha(mndwi.gt(-0.05)),
        "centro": area_ha(mndwi.gt(0.0)),
        "alto": area_ha(mndwi.gt(0.05)),
    }).getInfo()
    banda = (d["baixo"] - d["alto"]) / 2.0        # metade da amplitude = ± em torno do centro
    print("%4d| %9.1f | %6.1f | %10.1f | ±%6.1f" %
          (ano, d["baixo"], d["centro"], d["alto"], banda))
    linhas.append([ano, round(d["centro"], 1), round(banda, 1)])

with open("incerteza.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["ano", "area_ha", "incerteza_ha"])
    w.writerows(linhas)

media_banda = sum(l[2] for l in linhas) / len(linhas)
print("\nBanda de incerteza média: ±%.1f ha  (%.1f%% da área típica)"
      % (media_banda, 100 * media_banda / (sum(l[1] for l in linhas) / len(linhas))))
print("Guardado: incerteza.csv")
