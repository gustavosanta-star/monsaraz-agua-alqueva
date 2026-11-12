"""
Monsaraz em Mudança — passo 17
INDICADOR DE QUALIDADE DA IMAGEM QUE NÃO DEPENDE DA COTA ("teste do núcleo").

Problema: ~8 % das imagens de data única dão uma área errada (neblina, fumo, nuvens
finas, reflexo do sol) e, até aqui, só eram identificadas por comparação com a cota
oficial. Isso é circular e não serve para automatizar.

Ideia: o centro da albufeira é SEMPRE água. Se numa imagem "sem nuvens" esse núcleo
não aparecer como água, o problema é da imagem, não da albufeira.

Núcleo de água permanente (independente do Sentinel-2 e da cota):
  - JRC Global Surface Water, Yearly History: pixels classificados como "água
    permanente" em TODOS os anos de 2010 a 2021 (com pelo menos 8 anos observados);
  - encolhido 60 m para dentro, para ficar longe da margem;
  - dentro da envolvente da albufeira (passo 14).

Para cada data utilizável (as mesmas dos passos 12/16: lacuna < 0,2 % e cota medida):
  - nucleo_valido_pct  % do núcleo com dados válidos (máscara SCL);
  - nucleo_agua_pct    % do núcleo válido detetado como água (MNDWI > 0)  <- o indicador;
  - mndwi_nucleo       MNDWI médio no núcleo válido (medida contínua de "limpeza").

Uso:
  python 17_qualidade_imagem.py alqueva        uma albufeira
  python 17_qualidade_imagem.py                todas as que já têm dados/<id>/area_cota.csv
Se for interrompido, basta voltar a correr: continua de onde ficou.

Produz dados/<id>/qualidade.csv
"""

import csv
import json
import os
import sys
import time

import ee

ee.Initialize(project="projetobeta-513108")

PROJ = "EPSG:3763"
ESCALA = 20
LACUNA_FRAC = 50 / 24596
ENCOLHER_M = 60
LOTE = 12                       # datas por pedido
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()
MAU = [0, 1, 3, 8, 9, 10]

with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)
with open(os.path.join("dados", "resumo_envolventes.csv"), encoding="utf-8") as fh:
    ENV = {r["id"]: float(r["area_envolvente_ha"]) for r in csv.DictReader(fh)}

# --- núcleo de água permanente (JRC Yearly History, 2010–2021) ----------------
yh = ee.ImageCollection("JRC/GSW1_4/YearlyHistory").filterDate("2010-01-01", "2022-01-01")
if yh.size().getInfo() == 0:                      # recurso: filtrar pela propriedade "year"
    yh = (ee.ImageCollection("JRC/GSW1_4/YearlyHistory")
          .filter(ee.Filter.gte("year", 2010)).filter(ee.Filter.lte("year", 2021)))
n_anos = yh.size().getInfo()
if n_anos == 0:
    sys.exit("Não consegui filtrar o JRC Yearly History por ano — manda-me esta mensagem.")
print("JRC Yearly History: %d anos (2010–2021)" % n_anos)
yh = yh.select("waterClass")                      # 0 sem dados, 1 não-água, 2 sazonal, 3 permanente
permanente = yh.map(lambda im: im.eq(3)).sum()
observado = yh.map(lambda im: im.gte(1)).sum()
NUCLEO_BRUTO = permanente.eq(observado).And(observado.gte(8)).unmask(0)


def preparar(env):
    mascara_env = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(env)]), 1)
    nucleo = (NUCLEO_BRUTO.And(mascara_env)
              .focalMin(radius=ENCOLHER_M, units="meters", kernelType="circle").selfMask())
    caixa = env.bounds(1)
    ha = ee.Number(nucleo.multiply(pix).reduceRegion(
        ee.Reducer.sum(), caixa, ESCALA, PROJ, maxPixels=1e10, tileScale=4).values().get(0)).divide(1e4).getInfo()
    return nucleo, caixa, ha


def stats_data(data_str, nucleo, caixa):
    d = ee.Date(data_str)
    mos = (S2.filterBounds(caixa).filterDate(d, d.advance(1, "day"))
           .select(["B3", "B11", "SCL"]).mosaic())
    valido = mos.select("SCL").remap(MAU, [0] * len(MAU), 1)
    mndwi = mos.normalizedDifference(["B3", "B11"])
    agua = mndwi.gt(0).And(valido)
    img = (ee.Image.cat(valido.rename("v"), agua.rename("w"), mndwi.multiply(valido).rename("m"))
           .unmask(0).multiply(pix).updateMask(nucleo))
    r = img.reduceRegion(ee.Reducer.sum(), caixa, ESCALA, PROJ, maxPixels=1e10, tileScale=8)
    return ee.Feature(None, {"data": data_str, "v": r.get("v"), "w": r.get("w"), "m": r.get("m")})


def lote(datas, nucleo, caixa):
    fc = ee.FeatureCollection([stats_data(d, nucleo, caixa) for d in datas])
    return [f["properties"] for f in fc.getInfo()["features"]]


def lote_robusto(datas, nucleo, caixa):
    try:
        return lote(datas, nucleo, caixa)
    except Exception as e:
        print("    lote pesado (%s) — data a data" % str(e)[:40], flush=True)
    res = []
    for d in datas:
        for k in range(3):
            try:
                res.append(stats_data(d, nucleo, caixa).getInfo()["properties"])
                break
            except Exception:
                time.sleep(5)
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
        datas = [r["data"] for r in csv.DictReader(fh)
                 if r["cota_m"] and float(r["lacuna_ha"]) < LACUNA_FRAC * env_ha]

    with open(os.path.join(pasta, "envolvente.geojson"), encoding="utf-8") as fh:
        env = ee.Geometry(json.load(fh)["features"][0]["geometry"], None, False)
    nucleo, caixa, nucleo_ha = preparar(env)
    print("  núcleo de água permanente: %.0f ha (%.0f%% da envolvente) | datas a avaliar: %d"
          % (nucleo_ha, 100 * nucleo_ha / env_ha, len(datas)), flush=True)
    if nucleo_ha <= 0:
        print("  !! núcleo vazio — albufeira demasiado estreita para este teste")
        continue

    destino = os.path.join(pasta, "qualidade.csv")
    feitas = set()
    if os.path.exists(destino):
        with open(destino, encoding="utf-8") as fh:
            feitas = {r["data"] for r in csv.DictReader(fh)}
        print("  a retomar: %d datas já feitas" % len(feitas))
    else:
        with open(destino, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(["data", "nucleo_ha", "nucleo_valido_pct", "nucleo_agua_pct", "mndwi_nucleo"])

    por_fazer = [d for d in datas if d not in feitas]
    for i in range(0, len(por_fazer), LOTE):
        res = lote_robusto(por_fazer[i:i + LOTE], nucleo, caixa)
        with open(destino, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            for r in res:
                v, a, m = (r.get("v") or 0) / 1e4, (r.get("w") or 0) / 1e4, (r.get("m") or 0) / 1e4
                w.writerow([r["data"], round(nucleo_ha, 1), round(100 * v / nucleo_ha, 2),
                            round(100 * a / v, 2) if v > 0 else "", round(m / v, 4) if v > 0 else ""])
        print("  %d / %d datas" % (min(i + LOTE, len(por_fazer)), len(por_fazer)), flush=True)

    with open(destino, encoding="utf-8") as fh:
        q = [float(r["nucleo_agua_pct"]) for r in csv.DictReader(fh) if r["nucleo_agua_pct"]]
    if q:
        q.sort()
        print("  % do núcleo visto como água: mediana {:.1f} | mínimo {:.1f} | datas abaixo de 98 %: {}".format(
            q[len(q) // 2], q[0], sum(1 for x in q if x < 98)))
    print("  guardado: %s" % destino)
