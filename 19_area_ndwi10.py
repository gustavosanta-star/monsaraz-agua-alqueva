"""
Monsaraz em Mudança — passo 19 (Fase 3)
ÁREA DE ÁGUA A 10 m (NDWI) NAS DATAS UTILIZÁVEIS DE CADA ALBUFEIRA.

O passo 18 mostrou que o método original (MNDWI a 20 m, limiar 0) perde a faixa de margem
(águas rasas, turvas, com vegetação): no Alqueva mede 90 % da área implícita no volume
oficial; nas pequenas, 53–76 %. O NDWI = (B3 − B8)/(B3 + B8) usa duas bandas nativas
de 10 m e, com limiar −0,1, deu 1,00 no Alqueva e 0,98 no Caia.

O limiar −0,1 foi escolhido SÓ com o Alqueva e o Caia (as duas albufeiras cuja curva oficial
é coerente com o máximo histórico do JRC) e fica FIXO para as outras: assim as outras são um
teste e não uma calibração.

Mede, nas mesmas datas usadas nos passos 16/17 (lacuna < 0,2 % da envolvente, cota medida):
  agua_ha    área NDWI > −0,1 com dados válidos (SCL), a 10 m, dentro da envolvente
  valido_ha  área com dados válidos
Junta-lhe o indicador do núcleo do passo 17 (nucleo_agua_pct).

Uso:
  python 19_area_ndwi10.py               todas
  python 19_area_ndwi10.py caia enxoe    só estas
Se for interrompido, basta voltar a correr: continua de onde ficou.

Produz dados/<id>/area_cota_ndwi10.csv com as colunas do area_cota.csv + nucleo_agua_pct
"""

import csv
import json
import os
import sys
import time

import ee

ee.Initialize(project="projetobeta-513108")

PROJ = "EPSG:3763"
LIMIAR = -0.1
LACUNA_FRAC = 50 / 24596
LOTE = {"alqueva": 4}           # datas por pedido (as outras: 12)
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()
MAU = [0, 1, 3, 8, 9, 10]

with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)
with open(os.path.join("dados", "resumo_envolventes.csv"), encoding="utf-8") as fh:
    ENV = {r["id"]: float(r["area_envolvente_ha"]) for r in csv.DictReader(fh)}


def stats_data(data_str, mascara, caixa, ts):
    d = ee.Date(data_str)
    mos = (S2.filterBounds(caixa).filterDate(d, d.advance(1, "day"))
           .select(["B3", "B8", "SCL"]).mosaic())
    valido = mos.select("SCL").remap(MAU, [0] * len(MAU), 1)
    agua = mos.normalizedDifference(["B3", "B8"]).gt(LIMIAR).And(valido)
    img = (ee.Image.cat(valido.rename("v"), agua.rename("w")).unmask(0)
           .multiply(pix).updateMask(mascara))
    r = img.reduceRegion(ee.Reducer.sum(), caixa, 10, PROJ, maxPixels=1e11, tileScale=ts)
    return ee.Feature(None, {"data": data_str, "v": r.get("v"), "w": r.get("w")})


def lote_robusto(datas, mascara, caixa, ts):
    try:
        fc = ee.FeatureCollection([stats_data(d, mascara, caixa, ts) for d in datas])
        return [f["properties"] for f in fc.getInfo()["features"]]
    except Exception as e:
        print("    lote pesado (%s) — data a data" % str(e)[:40], flush=True)
    res = []
    for d in datas:
        for k in range(3):
            try:
                res.append(stats_data(d, mascara, caixa, 16).getInfo()["properties"])
                break
            except Exception:
                time.sleep(8)
    return res


pedidas = sys.argv[1:]
for alb in CFG["albufeiras"]:
    pasta = os.path.join("dados", alb["id"])
    origem = os.path.join(pasta, "area_cota.csv")
    if (pedidas and alb["id"] not in pedidas) or not os.path.exists(origem):
        continue
    print("\n%s" % alb["nome"], flush=True)
    env_ha = ENV[alb["id"]]
    with open(origem, encoding="utf-8") as fh:
        base = [r for r in csv.DictReader(fh)
                if r["cota_m"] and float(r["lacuna_ha"]) < LACUNA_FRAC * env_ha]
    nucleo = {}
    q = os.path.join(pasta, "qualidade.csv")
    if os.path.exists(q):
        with open(q, encoding="utf-8") as fh:
            nucleo = {r["data"]: r["nucleo_agua_pct"] for r in csv.DictReader(fh)}

    with open(os.path.join(pasta, "envolvente.geojson"), encoding="utf-8") as fh:
        env = ee.Geometry(json.load(fh)["features"][0]["geometry"], None, False)
    mascara = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(env)]), 1).selfMask()
    caixa = env.bounds(1)
    ts = 16 if alb["id"] == "alqueva" else 4

    destino = os.path.join(pasta, "area_cota_ndwi10.csv")
    feitas = set()
    if os.path.exists(destino):
        with open(destino, encoding="utf-8") as fh:
            feitas = {r["data"] for r in csv.DictReader(fh)}
        print("  a retomar: %d datas já feitas" % len(feitas))
    else:
        with open(destino, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(["data", "agua_ha", "valido_ha", "lacuna_ha", "lacuna_pct",
                                     "n_granulos", "cota_m", "dif_dias", "agua_mndwi20_ha", "nucleo_agua_pct"])
    por_fazer = [r for r in base if r["data"] not in feitas]
    print("  datas a medir: %d" % len(por_fazer), flush=True)
    n = LOTE.get(alb["id"], 12)
    for i in range(0, len(por_fazer), n):
        grupo = por_fazer[i:i + n]
        res = {r["data"]: r for r in lote_robusto([g["data"] for g in grupo], mascara, caixa, ts)}
        with open(destino, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            for g in grupo:
                r = res.get(g["data"])
                if not r:
                    continue
                v, a = (r.get("v") or 0) / 1e4, (r.get("w") or 0) / 1e4
                w.writerow([g["data"], round(a, 2), round(v, 2), g["lacuna_ha"], g["lacuna_pct"],
                            g["n_granulos"], g["cota_m"], g["dif_dias"], g["agua_ha"],
                            nucleo.get(g["data"], "")])
        print("  %d / %d" % (min(i + n, len(por_fazer)), len(por_fazer)), flush=True)
    print("  guardado: %s" % destino)

print("\nFeito. Diz ao Claude para analisar (passo 20).")
