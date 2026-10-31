"""
Monsaraz em Mudança — passo 09
SEPARAR O ALQUEVA DE OUTRAS MASSAS DE ÁGUA (charcas, pequenas barragens).

Problema: o MNDWI deteta toda a água, incluindo charcas e pequenas barragens
agrícolas dentro da freguesia. O indicador é a água DO ALQUEVA, por isso essas
manchas não devem contar.

Método — envolvente da albufeira:
  1. JRC Global Surface Water, banda "max_extent": pixéis que estiveram com água
     pelo menos uma vez entre 1984 e 2021 (inclui o Alqueva no nível máximo).
  2. Numa região alargada (freguesia + 5 km), ficar só com as manchas contínuas
     grandes (>= 1024 pixéis de 30 m, ~92 ha). O Alqueva é uma mancha enorme e
     contínua; as charcas são pequenas e isoladas.
  3. Alargar a envolvente 1 pixel (30 m) por segurança.
  4. Água do Alqueva = água detetada (MNDWI > 0) DENTRO da envolvente.
     Uma poça isolada no leito do Alqueva, a níveis baixos, continua a contar.

Este script só MEDE o efeito; ainda não altera os resultados anteriores.
Produz:
  - agua_fora_alqueva.csv  (por ano: total, Alqueva, excluída)
  - envolvente_2025.png    (azul = Alqueva; magenta = excluída)
"""

import ee
import json
import csv
import urllib.request

ee.Initialize(project="projetobeta-513108")

with open("monsaraz.geojson", encoding="utf-8") as fh:
    geom = json.load(fh)["features"][0]["geometry"]
aoi = ee.Geometry(geom, None, False)
regiao = aoi.buffer(5000)

PROJ = "EPSG:3763"
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()

# --- envolvente da albufeira (a partir da extensão máxima do JRC) ----------
jrc = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("max_extent")
proj30 = jrc.projection()
maxext = jrc.unmask(0).clip(regiao).reproject(proj30)
tamanho = (maxext.selfMask()
           .connectedPixelCount(maxSize=1024, eightConnected=True)
           .reproject(proj30))
envolvente = (tamanho.gte(1024).unmask(0)
              .focal_max(1, "square", "pixels").reproject(proj30)
              .rename("env"))


def mascara(img):
    scl = img.select("SCL")
    mau = (scl.eq(0).Or(scl.eq(1)).Or(scl.eq(3)).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)))
    return img.updateMask(mau.Not()).select(["B2", "B3", "B4", "B11"])


def comp_t3(ano):
    return (S2.filterBounds(aoi).filterDate(f"{ano}-07-01", f"{ano}-09-30")
            .map(mascara).median())


def area_ha(mask):
    m2 = (mask.selfMask().multiply(pix)
          .reduceRegion(ee.Reducer.sum(), aoi, 10, PROJ, maxPixels=1e9).values().get(0))
    return ee.Number(ee.Algorithms.If(m2, m2, 0)).divide(10000)


linhas = []
print("ano | total (ha) | Alqueva (ha) | excluída (ha)")
print("----+------------+--------------+--------------")
for ano in range(2017, 2026):
    agua = comp_t3(ano).normalizedDifference(["B3", "B11"]).gt(0).unmask(0).clip(aoi)
    d = ee.Dictionary({
        "total": area_ha(agua),
        "alq": area_ha(agua.And(envolvente)),
        "fora": area_ha(agua.And(envolvente.Not())),
    }).getInfo()
    print("%4d| %10.1f | %12.1f | %12.1f" % (ano, d["total"], d["alq"], d["fora"]))
    linhas.append([ano, round(d["total"], 1), round(d["alq"], 1), round(d["fora"], 1)])

with open("agua_fora_alqueva.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["ano", "total_ha", "alqueva_ha", "excluida_ha"])
    w.writerows(linhas)

# --- mapa de verificação (2025) ---------------------------------------------
comp = comp_t3(2025).clip(aoi)
agua = comp.normalizedDifference(["B3", "B11"]).gt(0)
rgb = comp.visualize(bands=["B4", "B3", "B2"], min=0, max=3000)
vis = (rgb
       .blend(agua.And(envolvente).selfMask().visualize(palette=["1f77b4"]))
       .blend(agua.And(envolvente.Not()).selfMask().visualize(palette=["ff00ff"])))
cont = ee.Image().paint(ee.FeatureCollection([ee.Feature(aoi)]), 0, 2)
vis = vis.blend(cont.visualize(palette=["FFFF00"]))
url = vis.getThumbURL({"region": aoi.bounds(), "dimensions": 1600, "format": "png"})
urllib.request.urlretrieve(url, "envolvente_2025.png")

print("\nGuardado: agua_fora_alqueva.csv")
print("Guardado: envolvente_2025.png (azul = Alqueva, magenta = excluída)")
