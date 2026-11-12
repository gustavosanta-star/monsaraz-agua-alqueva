"""
Monsaraz em Mudança — passo 15 (Fase 3)
ÁREA DE ÁGUA POR DATA EM CADA ALBUFEIRA + COTA OFICIAL DO MESMO DIA.

É o passo 11 generalizado (o mesmo método, sem alterações):
  - Sentinel-2 L2A, imagens de DATAS ÚNICAS (mosaico do dia), máscara de nuvens SCL;
  - água = MNDWI > 0, contada dentro da envolvente do passo 14;
  - "lacuna" = área da envolvente sem dados válidos (nuvens, sombras, falhas);
  - junta-se a cota oficial (SNIRH, passo 13) do mesmo dia, até 3 dias de diferença.

Para o Alqueva reaproveita o area_cota_alqueva.csv do passo 11 (não volta a calcular).

Corre ano a ano por albufeira (se um ano falhar, mês a mês; depois data a data) e guarda
o progresso em dados/<id>/area_parcial.csv: se for interrompido, basta voltar a correr.

Uso:
  python 15_area_cota_multi.py            todas
  python 15_area_cota_multi.py caia vigia só estas

Produz dados/<id>/area_cota.csv:
  data, agua_ha, valido_ha, lacuna_ha, lacuna_pct, n_granulos, cota_m, dif_dias
"""

import csv
import datetime as dt
import json
import os
import shutil
import sys
import time

import ee

ee.Initialize(project="projetobeta-513108")

PROJ = "EPSG:3763"
ESCALA = 20
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()
MAU = [0, 1, 3, 8, 9, 10]
MAX_DIF_DIAS = 3
ANOS = range(2017, 2026)

with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)
pedidas = sys.argv[1:]


def preparar(env):
    mascara = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(env)]), 1).selfMask()
    caixa = env.bounds(1)
    env_ha = ee.Number(mascara.multiply(pix).reduceRegion(
        ee.Reducer.sum(), caixa, ESCALA, PROJ, maxPixels=1e10).values().get(0)).divide(1e4).getInfo()
    return mascara, caixa, env_ha


def stats_data(data_str, mascara, caixa):
    d = ee.Date(data_str)
    imgs = S2.filterBounds(caixa).filterDate(d, d.advance(1, "day")).select(["B3", "B11", "SCL"])
    mos = imgs.mosaic()
    valido = mos.select("SCL").remap(MAU, [0] * len(MAU), 1)
    agua = mos.normalizedDifference(["B3", "B11"]).gt(0).And(valido)
    img = (ee.Image.cat(valido.rename("v"), agua.rename("w")).unmask(0)
           .multiply(pix).updateMask(mascara))
    r = img.reduceRegion(ee.Reducer.sum(), caixa, ESCALA, PROJ, maxPixels=1e10, tileScale=4)
    return ee.Feature(None, {"data": data_str,
                             "valido_ha": ee.Number(r.get("v")).divide(1e4),
                             "agua_ha": ee.Number(r.get("w")).divide(1e4),
                             "n_granulos": imgs.size()})


def datas_periodo(ini, fim, caixa):
    col = (S2.filterBounds(caixa).filterDate(ini, fim)
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40)))
    return (col.aggregate_array("system:time_start")
            .map(lambda t: ee.Date(t).format("YYYY-MM-dd")).distinct())


def processar(ini, fim, mascara, caixa):
    fc = ee.FeatureCollection(datas_periodo(ini, fim, caixa).map(lambda s: stats_data(s, mascara, caixa)))
    return [f["properties"] for f in fc.getInfo()["features"]]


def periodo_robusto(ano, mascara, caixa):
    """Um ano num pedido; se falhar, mês a mês; se um mês falhar, data a data."""
    try:
        return processar("%d-01-01" % ano, "%d-01-01" % (ano + 1), mascara, caixa)
    except Exception as e:
        print("    %d: ano pesado (%s) — mês a mês" % (ano, str(e)[:40]), flush=True)
    res = []
    for m in range(1, 13):
        ini = dt.date(ano, m, 1).isoformat()
        fim = dt.date(ano + (m == 12), m % 12 + 1, 1).isoformat()
        try:
            res += processar(ini, fim, mascara, caixa)
        except Exception:
            for data in datas_periodo(ini, fim, caixa).getInfo():
                for k in range(3):
                    try:
                        res.append(stats_data(data, mascara, caixa).getInfo()["properties"])
                        break
                    except Exception:
                        time.sleep(5)
    return res


def ler_cota(pasta):
    cota = {}
    with open(os.path.join(pasta, "snirh_diario.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["cota_m"]:
                cota[dt.date.fromisoformat(r["data"])] = float(r["cota_m"])
    return cota


def cota_mais_proxima(cota, d):
    for k in range(MAX_DIF_DIAS + 1):
        for s in ((0,) if k == 0 else (-k, k)):
            dd = d + dt.timedelta(days=s)
            if dd in cota:
                return cota[dd], s
    return None, None


for alb in CFG["albufeiras"]:
    if pedidas and alb["id"] not in pedidas:
        continue
    pasta = os.path.join("dados", alb["id"])
    saida_csv = os.path.join(pasta, "area_cota.csv")
    print("\n%s" % alb["nome"], flush=True)

    if alb["id"] == "alqueva" and os.path.exists("area_cota_alqueva.csv"):
        # passo 11: mesmas colunas, falta só a lacuna em %
        with open("area_cota_alqueva.csv", encoding="utf-8") as fh:
            linhas = list(csv.DictReader(fh))
        env_ha = float(linhas[0]["valido_ha"]) + float(linhas[0]["lacuna_ha"])
        with open(saida_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["data", "agua_ha", "valido_ha", "lacuna_ha", "lacuna_pct", "n_granulos", "cota_m", "dif_dias"])
            for r in linhas:
                w.writerow([r["data"], r["agua_ha"], r["valido_ha"], r["lacuna_ha"],
                            round(100 * float(r["lacuna_ha"]) / env_ha, 2), r["n_granulos"], r["cota_m"], r["dif_dias"]])
        print("  reaproveitado do passo 11: %d datas" % len(linhas))
        continue

    with open(os.path.join(pasta, "envolvente.geojson"), encoding="utf-8") as fh:
        env = ee.Geometry(json.load(fh)["features"][0]["geometry"], None, False)
    mascara, caixa, env_ha = preparar(env)
    print("  envolvente: %.0f ha" % env_ha, flush=True)

    parcial = os.path.join(pasta, "area_parcial.csv")
    feitos, linhas = set(), []
    if os.path.exists(parcial):
        with open(parcial, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                feitos.add(r["ano"])
                if r["data"]:
                    linhas.append({"data": r["data"], "agua_ha": float(r["agua_ha"]),
                                   "valido_ha": float(r["valido_ha"]), "n_granulos": int(r["n_granulos"])})
        print("  a retomar: %d anos já feitos" % len(feitos))
    else:
        with open(parcial, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(["ano", "data", "agua_ha", "valido_ha", "n_granulos"])

    for ano in ANOS:
        if str(ano) in feitos:
            continue
        res = periodo_robusto(ano, mascara, caixa)
        with open(parcial, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            if not res:
                w.writerow([ano, "", "", "", ""])
            for r in res:
                w.writerow([ano, r["data"], r["agua_ha"], r["valido_ha"], r["n_granulos"]])
        linhas += res
        limpas = sum(1 for r in res if env_ha - r["valido_ha"] <= 0.01 * env_ha)
        print("  %d: %3d datas, %3d com < 1%% de lacuna" % (ano, len(res), limpas), flush=True)

    cota = ler_cota(pasta)
    n_ok = 0
    with open(saida_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["data", "agua_ha", "valido_ha", "lacuna_ha", "lacuna_pct", "n_granulos", "cota_m", "dif_dias"])
        for r in sorted(linhas, key=lambda x: x["data"]):
            h, dif = cota_mais_proxima(cota, dt.date.fromisoformat(r["data"]))
            lac = max(env_ha - r["valido_ha"], 0)
            n_ok += (lac <= 0.01 * env_ha and h is not None)
            w.writerow([r["data"], round(r["agua_ha"], 2), round(r["valido_ha"], 2), round(lac, 2),
                        round(100 * lac / env_ha, 2), r["n_granulos"],
                        "" if h is None else round(h, 3), "" if dif is None else dif])
    print("  total: %d datas | com < 1%% de lacuna e cota: %d -> %s" % (len(linhas), n_ok, saida_csv), flush=True)

print("\nFeito. Próximo passo: 16 (curvas área–cota e volume de todas as albufeiras).")
