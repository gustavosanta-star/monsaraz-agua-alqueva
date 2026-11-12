"""
Monsaraz em Mudança — passo 17b
TESTE DO NÚCLEO, VERSÃO 2: núcleo definido pelo próprio Sentinel-2 no período de nível mais baixo.

Porquê: o núcleo do passo 17 vinha do JRC ("água permanente" 2010–2021). No Caia, no outono
de 2019, a albufeira desceu abaixo desse núcleo e o filtro rejeitou todas as datas de nível
mais baixo — rejeitava uma seca verdadeira, não imagens más.

Núcleo novo (sem JRC e sem cota):
  1. Na série de áreas do passo 19, procura-se a janela de 60 dias com a menor área mediana
     (mínimo 3 datas; a mediana aguenta uma imagem má isolada). É o período de nível mais baixo.
  2. Faz-se a imagem mediana do Sentinel-2 nessa janela (com máscara de nuvens).
  3. Núcleo = pixels que são água pelos DOIS índices (MNDWI > 0 e NDWI > 0), dentro da
     envolvente, encolhidos para dentro 60 m (ou 40 / 20 m se a albufeira for pequena demais).
  O contorno fica gravado em dados/<id>/nucleo.geojson e é reutilizado.

Limitação conhecida: se um dia a albufeira descer abaixo do mínimo de 2017–2025, o núcleo tem
de ser recalculado (o teste voltaria a confundir seca com imagem má).

Para cada data utilizável mede, dentro do núcleo:
  nucleo_agua_pct   % visto como água pelo MNDWI > 0 a 20 m   (o indicador do passo 17)
  nucleo_ndwi_pct   % visto como água pelo NDWI > −0,1 a 10 m  (o índice do passo 19)
  mndwi_nucleo      MNDWI médio

Requer o passo 19 (dados/<id>/area_cota_ndwi10.csv).
Uso:  python 17b_nucleo_s2.py            todas
      python 17b_nucleo_s2.py caia       só esta
Produz dados/<id>/nucleo.geojson e dados/<id>/qualidade_s2.csv
"""

import csv
import datetime as dt
import json
import os
import sys
import time
from statistics import median

import ee

ee.Initialize(project="projetobeta-513108")

PROJ = "EPSG:3763"
JANELA_DIAS = 60
MIN_DATAS = 3
EROSOES = [60, 40, 20]
NUCLEO_MIN_HA = 2
LOTE = {"alqueva": 6}
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
pix = ee.Image.pixelArea()
MAU = [0, 1, 3, 8, 9, 10]

with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)


def janela_minima(linhas):
    """Janela de JANELA_DIAS com a menor área mediana. linhas: [(data, área)] ordenadas."""
    melhor = None
    for i, (d0, _) in enumerate(linhas):
        fim = d0 + dt.timedelta(days=JANELA_DIAS)
        dentro = [a for d, a in linhas[i:] if d <= fim]
        if len(dentro) >= MIN_DATAS:
            m = median(dentro)
            if melhor is None or m < melhor[0]:
                melhor = (m, d0, fim, len(dentro))
    return melhor


def area_ha(img, caixa, escala=20):
    return ee.Number(img.multiply(pix).reduceRegion(
        ee.Reducer.sum(), caixa, escala, PROJ, maxPixels=1e10, tileScale=8).values().get(0)).divide(1e4).getInfo()


def construir_nucleo(alb, env, caixa, linhas):
    m, d0, d1, n = janela_minima(linhas)
    print("  período de nível mais baixo: %s a %s (%d datas, área mediana %.0f ha)" % (d0, d1, n, m), flush=True)
    comp = (S2.filterBounds(caixa).filterDate(d0.isoformat(), (d1 + dt.timedelta(days=1)).isoformat())
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
            .select(["B3", "B8", "B11", "SCL"])
            .map(lambda im: im.updateMask(im.select("SCL").remap(MAU, [0] * len(MAU), 1)))
            .median())
    agua = (comp.normalizedDifference(["B3", "B11"]).gt(0)
            .And(comp.normalizedDifference(["B3", "B8"]).gt(0)).unmask(0))
    dentro = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(env)]), 1)
    base = agua.And(dentro)
    for raio in EROSOES:
        nucleo = base.focalMin(radius=raio, units="meters", kernelType="circle").selfMask()
        ha = area_ha(nucleo, caixa)
        print("    encolhido %d m: %.1f ha" % (raio, ha), flush=True)
        if ha >= NUCLEO_MIN_HA:
            break
    else:
        return None, None, 0, (d0, d1)
    geo = None
    try:
        vec = nucleo.reduceToVectors(geometry=caixa, scale=20, crs=PROJ, geometryType="polygon",
                                     eightConnected=False, maxPixels=1e10, tileScale=8)
        geo = vec.geometry().simplify(10).getInfo()
    except Exception as e:
        print("    (não consegui gravar o contorno: %s — uso a imagem diretamente)" % str(e)[:60], flush=True)
    return nucleo, geo, raio, (d0, d1)


def stats_data(data_str, nucleo, caixa):
    d = ee.Date(data_str)
    mos = (S2.filterBounds(caixa).filterDate(d, d.advance(1, "day"))
           .select(["B3", "B8", "B11", "SCL"]).mosaic())
    valido = mos.select("SCL").remap(MAU, [0] * len(MAU), 1)
    mndwi = mos.normalizedDifference(["B3", "B11"])
    ndwi = mos.normalizedDifference(["B3", "B8"])
    img20 = (ee.Image.cat(valido.rename("v"), mndwi.gt(0).And(valido).rename("w"),
                          mndwi.multiply(valido).rename("m"))
             .unmask(0).multiply(pix).updateMask(nucleo))
    img10 = (ee.Image.cat(valido.rename("v10"), ndwi.gt(-0.1).And(valido).rename("n10"))
             .unmask(0).multiply(pix).updateMask(nucleo))
    r20 = img20.reduceRegion(ee.Reducer.sum(), caixa, 20, PROJ, maxPixels=1e10, tileScale=8)
    r10 = img10.reduceRegion(ee.Reducer.sum(), caixa, 10, PROJ, maxPixels=1e10, tileScale=8)
    return ee.Feature(None, r20.combine(r10).set("data", data_str))


def lote_robusto(datas, nucleo, caixa):
    try:
        fc = ee.FeatureCollection([stats_data(d, nucleo, caixa) for d in datas])
        return [f["properties"] for f in fc.getInfo()["features"]]
    except Exception as e:
        print("    lote pesado (%s) — data a data" % str(e)[:40], flush=True)
    res = []
    for d in datas:
        for k in range(3):
            try:
                res.append(stats_data(d, nucleo, caixa).getInfo()["properties"])
                break
            except Exception:
                time.sleep(8)
    return res


pedidas = sys.argv[1:]
for alb in CFG["albufeiras"]:
    pasta = os.path.join("dados", alb["id"])
    origem = os.path.join(pasta, "area_cota_ndwi10.csv")
    if (pedidas and alb["id"] not in pedidas) or not os.path.exists(origem):
        continue
    print("\n%s" % alb["nome"], flush=True)
    with open(origem, encoding="utf-8") as fh:
        base = list(csv.DictReader(fh))
    linhas = sorted((dt.date.fromisoformat(r["data"]), float(r["agua_ha"])) for r in base)
    with open(os.path.join(pasta, "envolvente.geojson"), encoding="utf-8") as fh:
        env = ee.Geometry(json.load(fh)["features"][0]["geometry"], None, False)
    caixa = env.bounds(1)

    f_nucleo = os.path.join(pasta, "nucleo.geojson")
    if os.path.exists(f_nucleo):
        with open(f_nucleo, encoding="utf-8") as fh:
            gj = json.load(fh)
        geo, raio = gj["features"][0]["geometry"], gj["features"][0]["properties"]["erosao_m"]
        nucleo = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(ee.Geometry(geo, None, False))]), 1).selfMask()
        print("  núcleo lido de nucleo.geojson")
    else:
        nucleo, geo, raio, (d0, d1) = construir_nucleo(alb, env, caixa, linhas)
        if nucleo is None:
            print("  !! albufeira pequena demais no nível mínimo: sem núcleo, sem teste")
            continue
        if geo is not None:
            with open(f_nucleo, "w", encoding="utf-8") as fh:
                json.dump({"type": "FeatureCollection", "features": [{
                    "type": "Feature",
                    "properties": {"nome": "Núcleo de água de %s" % alb["nome"], "erosao_m": raio,
                                   "periodo": "%s a %s" % (d0, d1)},
                    "geometry": geo}]}, fh)
            nucleo = ee.Image(0).paint(ee.FeatureCollection([ee.Feature(ee.Geometry(geo, None, False))]), 1).selfMask()
    nucleo_ha = area_ha(nucleo, caixa)
    print("  núcleo: %.1f ha (encolhido %d m) | datas a avaliar: %d" % (nucleo_ha, raio, len(base)), flush=True)

    destino = os.path.join(pasta, "qualidade_s2.csv")
    feitas = set()
    if os.path.exists(destino):
        with open(destino, encoding="utf-8") as fh:
            feitas = {r["data"] for r in csv.DictReader(fh)}
        print("  a retomar: %d datas já feitas" % len(feitas))
    else:
        with open(destino, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(["data", "nucleo_ha", "erosao_m", "nucleo_valido_pct",
                                     "nucleo_agua_pct", "nucleo_ndwi_pct", "mndwi_nucleo"])
    por_fazer = [r["data"] for r in base if r["data"] not in feitas]
    n = LOTE.get(alb["id"], 12)
    for i in range(0, len(por_fazer), n):
        res = lote_robusto(por_fazer[i:i + n], nucleo, caixa)
        with open(destino, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            for r in res:
                v, a, m = (r.get("v") or 0) / 1e4, (r.get("w") or 0) / 1e4, (r.get("m") or 0) / 1e4
                v10, n10 = (r.get("v10") or 0) / 1e4, (r.get("n10") or 0) / 1e4
                w.writerow([r["data"], round(nucleo_ha, 1), raio, round(100 * v / nucleo_ha, 2),
                            round(100 * a / v, 2) if v > 0 else "", round(100 * n10 / v10, 2) if v10 > 0 else "",
                            round(m / v, 4) if v > 0 else ""])
        print("  %d / %d" % (min(i + n, len(por_fazer)), len(por_fazer)), flush=True)
    with open(destino, encoding="utf-8") as fh:
        q = sorted(float(r["nucleo_agua_pct"]) for r in csv.DictReader(fh) if r["nucleo_agua_pct"])
    if q:
        print("  % do núcleo visto como água: mediana {:.1f} | datas abaixo de 98 %: {} de {}".format(
            q[len(q) // 2], sum(1 for x in q if x < 98), len(q)))

print("\nFeito. Diz ao Claude para analisar.")
