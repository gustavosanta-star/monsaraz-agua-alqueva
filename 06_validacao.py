"""
Monsaraz em Mudança — passo 06 (v2)
VALIDAÇÃO contra o JRC Global Surface Water (Pekel et al., 2016).

Correção face à v1: o JRC (derivado do Landsat) tem muitos meses com pixéis
"sem dados", sobretudo no inverno. Comparar nesses meses contaminava a métrica.
Agora:
  - usa-se só os meses em que o JRC observou > 80% da freguesia;
  - as duas águas são comparadas SOBRE OS MESMOS pixéis observados pelo JRC.
Tudo a 30 m (resolução do JRC).
"""

import ee
import json
import csv
import urllib.request

ee.Initialize(project="projetobeta-513108")

with open("monsaraz.geojson", encoding="utf-8") as fh:
    geom = json.load(fh)["features"][0]["geometry"]
aoi = ee.Geometry(geom, None, False)

PROJ = "EPSG:3763"
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
JRC = ee.ImageCollection("JRC/GSW1_4/MonthlyHistory")
pix = ee.Image.pixelArea()
aoi_ha = aoi.area(1, ee.Projection(PROJ)).divide(10000).getInfo()
COBERTURA_MIN = 0.80


def mascara(img):
    scl = img.select("SCL")
    mau = (scl.eq(0).Or(scl.eq(1)).Or(scl.eq(3)).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)))
    return img.updateMask(mau.Not()).select(["B3", "B11"])


def area_ha(mask):
    m2 = (mask.selfMask().multiply(pix)
          .reduceRegion(ee.Reducer.sum(), aoi, 30, PROJ, maxPixels=1e9).values().get(0))
    return ee.Number(ee.Algorithms.If(m2, m2, 0)).divide(10000)


def nossa_agua(ini, fim):
    comp = S2.filterBounds(aoi).filterDate(ini, fim).map(mascara).median()
    return comp.normalizedDifference(["B3", "B11"]).gt(0).unmask(0).clip(aoi)


# ---------- 1. Correlação temporal, só em meses JRC válidos ----------
pares = []
print("ano-mês | cobertura JRC | nossa (ha) | JRC (ha)")
print("--------+---------------+------------+---------")
for ano in range(2017, 2022):
    for m in range(1, 13):
        ini = ee.Date.fromYMD(ano, m, 1); fim = ini.advance(1, "month")
        jrc_img = JRC.filterDate(ini, fim).first()
        valido = ee.Image(jrc_img).select("water").neq(0)       # 1 = observado
        nossa = nossa_agua(ini, fim).updateMask(valido)
        jrc = ee.Image(jrc_img).select("water").eq(2).updateMask(valido)
        d = ee.Dictionary({"cov": area_ha(valido), "na": area_ha(nossa),
                           "ja": area_ha(jrc)}).getInfo()
        cob = d["cov"] / aoi_ha
        if cob >= COBERTURA_MIN:
            pares.append((ano, m, d["na"], d["ja"]))
            print("%4d-%02d |     %5.1f%%    | %10.1f | %8.1f" % (ano, m, 100*cob, d["na"], d["ja"]))

def pearson(x, y):
    n = len(x); mx = sum(x)/n; my = sum(y)/n
    sxy = sum((a-mx)*(b-my) for a, b in zip(x, y))
    sxx = sum((a-mx)**2 for a in x); syy = sum((b-my)**2 for b in y)
    return sxy/((sxx*syy)**0.5) if sxx and syy else float("nan")

xs = [p[2] for p in pares]; ys = [p[3] for p in pares]
r = pearson(xs, ys)
print("\nMeses JRC válidos usados: %d" % len(pares))
print("Correlação (Pearson r) nossa vs JRC: %.3f" % r)

with open("validacao.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["ano", "mes", "area_nossa_ha", "area_jrc_ha"])
    w.writerows(pares)

# ---------- 2. Concordância espacial (IoU) — julho 2020 ----------
ini, fim = ee.Date.fromYMD(2020, 7, 1), ee.Date.fromYMD(2020, 8, 1)
jrc_img = ee.Image(JRC.filterDate(ini, fim).first())
valido = jrc_img.select("water").neq(0)
nossa = nossa_agua(ini, fim).updateMask(valido)
jrc = jrc_img.select("water").eq(2).updateMask(valido)

inter = nossa.And(jrc); uniao = nossa.Or(jrc)
ai = area_ha(inter).getInfo(); au = area_ha(uniao).getInfo()
iou = ai/au if au else 0
print("\nConcordância espacial (julho 2020): IoU = %.1f%%" % (100*iou))

comp = (S2.filterBounds(aoi).filterDate(ini, fim)
        .map(lambda im: im.updateMask(im.select("SCL").remap([0,1,3,8,9,10],[0,0,0,0,0,0],1))).median())
rgb = comp.visualize(bands=["B4","B3","B2"], min=0, max=3000)
vis = (rgb
       .blend(inter.selfMask().visualize(palette=["1f77b4"]))
       .blend(nossa.And(jrc.Not()).selfMask().visualize(palette=["f5a623"]))
       .blend(jrc.And(nossa.Not()).selfMask().visualize(palette=["e040fb"])))
cont = ee.Image().paint(ee.FeatureCollection([ee.Feature(aoi)]), 0, 2)
vis = vis.blend(cont.visualize(palette=["FFFF00"]))
url = vis.getThumbURL({"region": aoi.bounds(), "dimensions": 1024, "format": "png"})
urllib.request.urlretrieve(url, "validacao_iou.png")
print("Guardado: validacao_iou.png (azul=ambos, laranja=só nossa, magenta=só JRC)")
print("Guardado: validacao.csv")
