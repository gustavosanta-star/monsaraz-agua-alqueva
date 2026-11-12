"""
Monsaraz em Mudança — passo 18 (Fase 3)
DE ONDE VEM A DIFERENÇA PARA A ÁREA OFICIAL? Teste de resolução e de limiar.

O passo 16 mostrou que a área do satélite fica abaixo da área implícita no volume oficial:
~8 % no Alqueva, mas 25–45 % nas albufeiras pequenas. Há duas explicações possíveis:
  (a) o satélite perde água na margem (pixels mistos de 20 m, braços estreitos, vegetação);
  (b) as curvas oficiais de volume das albufeiras pequenas sobrestimam a área
      (tabelas de projeto antigas, assoreamento).

Este passo separa as duas. Em ~10 datas limpas por albufeira (espalhadas pela gama de cotas),
mede a área de água de várias maneiras, dentro da mesma envolvente:
  - MNDWI a 20 m com limiares 0 (o método atual), −0,1, −0,2 e −0,3;
  - classe "água" da máscara SCL do Sentinel-2 (20 m);
  - NDWI = (B3 − B8)/(B3 + B8) a 10 m, limiares 0 e −0,1  (bandas nativas de 10 m);
  - MNDWI a 10 m com a B11 interpolada (bilinear), limiar 0.
E, uma vez por albufeira, a extensão máxima de água do JRC (Landsat, 1984–2021) dentro da envolvente.

Leitura: se as medidas a 10 m (ou limiares mais baixos) aproximarem a área da oficial -> (a).
Se todas as medidas do satélite concordarem entre si e com o máximo do JRC, e a oficial
continuar acima -> (b).

Requer o passo 16 já corrido (usa dados/<id>/area_cota_ajuste.csv para escolher datas limpas).
Uso:  python 18_resolucao.py            todas
      python 18_resolucao.py enxoe caia só estas
Produz dados/<id>/resolucao.csv e dados/resumo_jrc_max.csv
"""

import csv
import json
import os
import sys
import time

import ee

ee.Initialize(project="projetobeta-513108")

PROJ = "EPSG:3763"
N_DATAS = {"alqueva": 6}          # as outras: 10
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()
MAU = [0, 1, 3, 8, 9, 10]
JRC = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("max_extent")

with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)

COLS20 = ["mndwi20_0", "mndwi20_m01", "mndwi20_m02", "mndwi20_m03", "scl_agua"]
COLS10 = ["ndwi10_0", "ndwi10_m01", "mndwi10_0"]


def medir(data_str, mascara, caixa):
    d = ee.Date(data_str)
    imgs = S2.filterBounds(caixa).filterDate(d, d.advance(1, "day"))
    mos = imgs.select(["B3", "B8", "B11", "SCL"]).mosaic()
    b11_suave = imgs.map(lambda im: im.select("B11").resample("bilinear")).mosaic()
    valido = mos.select("SCL").remap(MAU, [0] * len(MAU), 1)
    mndwi = mos.normalizedDifference(["B3", "B11"])
    ndwi = mos.normalizedDifference(["B3", "B8"])
    mndwi10 = mos.select("B3").subtract(b11_suave).divide(mos.select("B3").add(b11_suave))

    img20 = ee.Image.cat(
        mndwi.gt(0).rename("mndwi20_0"), mndwi.gt(-0.1).rename("mndwi20_m01"),
        mndwi.gt(-0.2).rename("mndwi20_m02"), mndwi.gt(-0.3).rename("mndwi20_m03"),
        mos.select("SCL").eq(6).rename("scl_agua"),
    ).And(valido).unmask(0).multiply(pix).updateMask(mascara)
    img10 = ee.Image.cat(
        ndwi.gt(0).rename("ndwi10_0"), ndwi.gt(-0.1).rename("ndwi10_m01"),
        mndwi10.gt(0).rename("mndwi10_0"),
    ).And(valido).unmask(0).multiply(pix).updateMask(mascara)

    r20 = img20.reduceRegion(ee.Reducer.sum(), caixa, 20, PROJ, maxPixels=1e10, tileScale=4)
    r10 = img10.reduceRegion(ee.Reducer.sum(), caixa, 10, PROJ, maxPixels=1e10, tileScale=8)
    return r20.combine(r10)


pedidas = sys.argv[1:]
jrc_linhas = []
for alb in CFG["albufeiras"]:
    pasta = os.path.join("dados", alb["id"])
    ajuste = os.path.join(pasta, "area_cota_ajuste.csv")
    if (pedidas and alb["id"] not in pedidas) or not os.path.exists(ajuste):
        continue
    print("\n%s" % alb["nome"], flush=True)
    with open(os.path.join(pasta, "envolvente.geojson"), encoding="utf-8") as fh:
        env = ee.Geometry(json.load(fh)["features"][0]["geometry"], None, False)
    mascara = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(env)]), 1).selfMask()
    caixa = env.bounds(1)

    jrc_ha = ee.Number(JRC.selfMask().multiply(pix).updateMask(mascara).reduceRegion(
        ee.Reducer.sum(), caixa, 30, PROJ, maxPixels=1e10, tileScale=4).values().get(0)).divide(1e4).getInfo()
    jrc_linhas.append([alb["id"], alb["nome"], round(jrc_ha, 1)])
    print("  extensão máxima de água do JRC (1984–2021): %.0f ha" % jrc_ha, flush=True)

    with open(ajuste, encoding="utf-8") as fh:
        limpas = sorted((r for r in csv.DictReader(fh) if r["anomala"] == "nao"),
                        key=lambda r: float(r["cota_m"]))
    n = min(N_DATAS.get(alb["id"], 10), len(limpas))
    escolhidas = [limpas[round(i * (len(limpas) - 1) / (n - 1))] for i in range(n)] if n > 1 else limpas

    destino = os.path.join(pasta, "resolucao.csv")
    feitas = set()
    if os.path.exists(destino):
        with open(destino, encoding="utf-8") as fh:
            feitas = {r["data"] for r in csv.DictReader(fh)}
    else:
        with open(destino, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(["data", "cota_m", "passo15_ha"] + COLS20 + COLS10)

    for r in escolhidas:
        if r["data"] in feitas:
            continue
        res = None
        for k in range(3):
            try:
                res = medir(r["data"], mascara, caixa).getInfo()
                break
            except Exception as e:
                print("    %s: nova tentativa (%s)" % (r["data"], str(e)[:50]), flush=True)
                time.sleep(8)
        if res is None:
            print("    %s: falhou, a saltar" % r["data"])
            continue
        v = [round((res.get(c) or 0) / 1e4, 2) for c in COLS20 + COLS10]
        with open(destino, "a", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow([r["data"], r["cota_m"], r["agua_ha"]] + v)
        print("  %s cota %7.2f | MNDWI 20 m: %7.1f | limiar −0,2: %7.1f | NDWI 10 m: %7.1f | MNDWI 10 m: %7.1f ha"
              % (r["data"], float(r["cota_m"]), v[0], v[2], v[5], v[7]), flush=True)

if jrc_linhas:
    caminho = os.path.join("dados", "resumo_jrc_max.csv")
    antigos = {}
    if os.path.exists(caminho):
        with open(caminho, encoding="utf-8") as fh:
            antigos = {r[0]: r for r in list(csv.reader(fh))[1:]}
    for r in jrc_linhas:
        antigos[r[0]] = r
    with open(caminho, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "nome", "jrc_max_ha"])
        w.writerows(antigos.values())
print("\nGuardado: dados/<albufeira>/resolucao.csv e dados/resumo_jrc_max.csv")
