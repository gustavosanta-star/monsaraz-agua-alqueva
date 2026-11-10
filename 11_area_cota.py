"""
Monsaraz em Mudança — passo 11 (Fase 1)
CURVA ÁREA–COTA DO ALQUEVA: validação independente com a cota oficial (SNIRH).

Ideia: numa albufeira, a área inundada é uma função crescente do nível da água.
Se a área medida por satélite seguir de perto a cota medida no terreno, o método
de deteção da água está validado por uma fonte que não tem nada a ver com imagens.
É também a base da Fase 2 (volume: dV = A · dh).

Diferenças face aos passos anteriores:
  - Área do ALQUEVA INTEIRO (incluindo a parte espanhola), não só da freguesia,
    porque é a área da albufeira toda que se relaciona com a cota.
  - Imagens de DATAS ÚNICAS (não medianas trimestrais), cada uma com a cota do
    mesmo dia (ou até 3 dias de diferença).

Envolvente do Alqueva:
  - extensão máxima de água do JRC Global Surface Water (1984–2021);
  - corte a sul, logo a jusante da barragem, para excluir o Pedrógão;
  - fica-se com a maior massa de água contínua (o Alqueva) + 30 m de margem.

Requer: snirh_alqueva_diario.csv (passo 10) na mesma pasta.

Produz:
  - area_cota_parcial.csv      progresso mês a mês (permite retomar se for interrompido)
  - area_cota_alqueva.csv       uma linha por data: área, lacuna por nuvens, cota
  - alqueva_envolvente.geojson  contorno usado (simplificado)
  - envolvente_alqueva.png      verificação visual (contorno a amarelo, água a azul)
"""

import ee
import csv
import json
import datetime as dt
import urllib.request
import os

ee.Initialize(project="projetobeta-513108")

PROJ = "EPSG:3763"
ESCALA = 20                      # m — resolução nativa da B11 (a banda que limita o MNDWI)
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()
MAU = [0, 1, 3, 8, 9, 10]        # SCL: sem dados, saturado, sombra de nuvem, nuvens, cirrus
MAX_DIF_DIAS = 3                 # diferença máxima entre a imagem e a medição de cota

# --- 1. Envolvente do Alqueva ----------------------------------------------
# Corte sul em 38.199 N: a barragem está a ~38.1975 N; o Pedrógão fica a jusante (a sul).
# O contorno é calculado UMA vez e gravado em alqueva_envolvente.geojson; nas execuções
# seguintes é lido do ficheiro (recalculá-lo em cada pedido esgotava a memória do Earth Engine).
if os.path.exists("alqueva_envolvente.geojson"):
    with open("alqueva_envolvente.geojson", encoding="utf-8") as fh:
        envolvente = ee.Geometry(json.load(fh)["features"][0]["geometry"], None, False)
    print("Envolvente lida de alqueva_envolvente.geojson")
else:
    regiao = ee.Geometry.Rectangle([-7.75, 38.199, -6.90, 38.85])
    jrc = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("max_extent")
    vec = jrc.selfMask().reduceToVectors(
        geometry=regiao, scale=30, geometryType="polygon", eightConnected=True,
        labelProperty="agua", maxPixels=1e10)
    vec = vec.map(lambda f: f.set("area_m2", f.geometry().area(30)))
    alqueva = ee.Feature(vec.sort("area_m2", False).first()).geometry()
    geo = alqueva.buffer(30, 10).simplify(30).getInfo()
    with open("alqueva_envolvente.geojson", "w", encoding="utf-8") as fh:
        json.dump({"type": "FeatureCollection", "features": [{
            "type": "Feature",
            "properties": {"nome": "Envolvente do Alqueva (JRC max_extent, a montante da barragem)"},
            "geometry": geo}]}, fh)
    envolvente = ee.Geometry(geo, None, False)

# Máscara raster da envolvente (muito mais leve do que reduzir sobre o polígono complexo)
mascara_env = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(envolvente)]), 1).selfMask()
caixa = envolvente.bounds(1)
env_ha = (mascara_env.multiply(pix).reduceRegion(ee.Reducer.sum(), caixa, ESCALA, PROJ,
          maxPixels=1e10).values().get(0))
env_ha = ee.Number(env_ha).divide(1e4).getInfo()
print("Envolvente do Alqueva: %.0f ha" % env_ha)


# --- 2. Área de água por data única -----------------------------------------
def stats_data(data_str):
    d = ee.Date(data_str)
    imgs = S2.filterBounds(caixa).filterDate(d, d.advance(1, "day")).select(["B3", "B11", "SCL"])  # só as bandas usadas: evita o erro de coleção não homogénea
    mos = imgs.mosaic()
    valido = mos.select("SCL").remap(MAU, [0] * len(MAU), 1)       # 1 = pixel utilizável
    agua = mos.normalizedDifference(["B3", "B11"]).gt(0).And(valido)
    img = (ee.Image.cat(valido.rename("v"), agua.rename("w")).unmask(0)
           .multiply(pix).updateMask(mascara_env))
    r = img.reduceRegion(ee.Reducer.sum(), caixa, ESCALA, PROJ,
                         maxPixels=1e10, tileScale=8)
    return ee.Feature(None, {
        "data": data_str,
        "valido_ha": ee.Number(r.get("v")).divide(1e4),
        "agua_ha": ee.Number(r.get("w")).divide(1e4),
        "n_granulos": imgs.size(),
    })


def datas_periodo(ini, fim):
    col = (S2.filterBounds(caixa).filterDate(ini, fim)
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40)))
    return (col.aggregate_array("system:time_start")
            .map(lambda t: ee.Date(t).format("YYYY-MM-dd")).distinct())


def processar(ini, fim):
    fc = ee.FeatureCollection(datas_periodo(ini, fim).map(stats_data))
    return [f["properties"] for f in fc.getInfo()["features"]]


PARCIAL = "area_cota_parcial.csv"
CAMPOS = ["data", "agua_ha", "valido_ha", "n_granulos"]
feitos = set()
linhas = []
if os.path.exists(PARCIAL):                      # retomar uma execução interrompida
    with open(PARCIAL, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            feitos.add(r["mes"])
            if r["data"]:
                linhas.append({"data": r["data"], "agua_ha": float(r["agua_ha"]),
                               "valido_ha": float(r["valido_ha"]), "n_granulos": int(r["n_granulos"])})
    print("A retomar: %d meses já processados" % len(feitos))
else:
    with open(PARCIAL, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerow(["mes"] + CAMPOS)

for ano in range(2017, 2026):
    for m in range(1, 13):
        mes = "%d-%02d" % (ano, m)
        if mes in feitos:
            continue
        ini = dt.date(ano, m, 1)
        fim = dt.date(ano + (m == 12), m % 12 + 1, 1)
        try:
            res = processar(ini.isoformat(), fim.isoformat())
        except Exception as e:                   # mês pesado demais: uma data de cada vez
            print("  %s: a processar data a data (%s)" % (mes, str(e)[:50]), flush=True)
            res = []
            for data in datas_periodo(ini.isoformat(), fim.isoformat()).getInfo():
                res.append(stats_data(data).getInfo()["properties"])
        with open(PARCIAL, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            if not res:
                w.writerow([mes, "", "", "", ""])          # mês sem imagens: marcar como feito
            for r in res:
                w.writerow([mes, r["data"], r["agua_ha"], r["valido_ha"], r["n_granulos"]])
        linhas += res
        bons = sum(1 for r in res if env_ha - r["valido_ha"] < 50)
        print("  %s: %2d datas, %2d sem lacunas relevantes" % (mes, len(res), bons), flush=True)

# --- 3. Juntar a cota oficial (SNIRH) ---------------------------------------
cota = {}
with open("snirh_alqueva_diario.csv", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        if r["cota_m"]:
            cota[dt.date.fromisoformat(r["data"])] = float(r["cota_m"])


def cota_mais_proxima(d):
    for k in range(MAX_DIF_DIAS + 1):
        for s in ((0,) if k == 0 else (-k, k)):
            dd = d + dt.timedelta(days=s)
            if dd in cota:
                return cota[dd], s
    return None, None


saida = []
for r in sorted(linhas, key=lambda x: x["data"]):
    d = dt.date.fromisoformat(r["data"])
    h, dif = cota_mais_proxima(d)
    saida.append([r["data"], round(r["agua_ha"], 1), round(r["valido_ha"], 1),
                  round(env_ha - r["valido_ha"], 1), r["n_granulos"],
                  "" if h is None else round(h, 3), "" if dif is None else dif])

with open("area_cota_alqueva.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["data", "agua_ha", "valido_ha", "lacuna_ha", "n_granulos", "cota_m", "dif_dias"])
    w.writerows(saida)

uteis = [s for s in saida if s[3] < 50 and s[5] != ""]
print("\nDatas processadas: %d | sem lacunas e com cota: %d" % (len(saida), len(uteis)))
if len(uteis) > 2:
    xs = [s[5] for s in uteis]; ys = [s[1] for s in uteis]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs); syy = sum((b - my) ** 2 for b in ys)
    print("Pré-visualização: r (área vs cota) = %.3f" % (sxy / (sxx * syy) ** 0.5))
print("Guardado: area_cota_alqueva.csv, alqueva_envolvente.geojson")

# --- 4. Imagem de verificação ------------------------------------------------
comp = (S2.filterBounds(caixa).filterDate("2025-07-01", "2025-09-30")
        .select(["B2", "B3", "B4", "B11", "SCL"])
        .map(lambda im: im.updateMask(im.select("SCL").remap(MAU, [0] * len(MAU), 1)))
        .median())
rgb = comp.visualize(bands=["B4", "B3", "B2"], min=0, max=3000)
agua = comp.normalizedDifference(["B3", "B11"]).gt(0).updateMask(mascara_env)
cont = ee.Image().paint(ee.FeatureCollection([ee.Feature(envolvente)]), 0, 2)
vis = (rgb.blend(agua.selfMask().visualize(palette=["1f77b4"]))
       .blend(cont.visualize(palette=["FFFF00"])))
url = vis.getThumbURL({"region": caixa.buffer(3000), "dimensions": 1400, "format": "png"})
urllib.request.urlretrieve(url, "envolvente_alqueva.png")
print("Guardado: envolvente_alqueva.png (contorno amarelo = envolvente; água azul = dentro dela)")
