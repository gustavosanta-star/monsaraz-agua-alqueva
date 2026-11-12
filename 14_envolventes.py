"""
Monsaraz em Mudança — passo 14 (Fase 3)
ENVOLVENTE DE CADA ALBUFEIRA (albufeiras.json), a partir do JRC Global Surface Water.

Como no passo 11 para o Alqueva:
  - extensão máxima de água do JRC (1984–2021) numa caixa à volta da barragem;
  - polígonos de água a menos de 1 km da estação do SNIRH (que fica na barragem);
  - escolhe-se o MAIOR deles (a albufeira) + 30 m de margem.
Os outros candidatos ficam registados e aparecem a vermelho na imagem de verificação:
se o escolhido estiver errado (ex.: o rio a jusante), corrige-se no albufeiras.json
com "caixa": [lon_min, lat_min, lon_max, lat_max] (zona de procura) ou
"excluir": [[lon_min, lat_min, lon_max, lat_max], ...] (zonas a retirar) e volta-se a
correr só essa albufeira.

Uso:
  python 14_envolventes.py              todas (salta as que já têm envolvente)
  python 14_envolventes.py caia --refazer

Produz, em dados/<id>/:
  - envolvente.geojson   contorno usado
  - envolvente.png       verificação: amarelo = escolhido, vermelho = outros candidatos,
                         azul = água no verão de 2025, ponto magenta = estação SNIRH
E dados/resumo_envolventes.csv (área e perímetro de cada envolvente).
"""

import csv
import json
import os
import shutil
import sys
import urllib.request

import ee

ee.Initialize(project="projetobeta-513108")

PROJ = "EPSG:3763"
RAIO_ESTACAO = 1000          # m — candidatos a menos disto da estação
MEIA_CAIXA = (0.15, 0.12)    # graus (lon, lat) à volta da estação, se não houver "caixa"
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
MAU = [0, 1, 3, 8, 9, 10]
JRC = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("max_extent")

with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)

args = [a for a in sys.argv[1:] if not a.startswith("--")]
refazer = "--refazer" in sys.argv


def guardar_geojson(caminho, geo, nome):
    with open(caminho, "w", encoding="utf-8") as fh:
        json.dump({"type": "FeatureCollection", "features": [{
            "type": "Feature", "properties": {"nome": nome}, "geometry": geo}]}, fh)


def imagem(alb, env, outros, ponto, caminho):
    caixa = env.bounds(1).buffer(1500, 1)
    comp = (S2.filterBounds(caixa).filterDate("2025-07-01", "2025-09-30")
            .select(["B2", "B3", "B4", "B11", "SCL"])
            .map(lambda im: im.updateMask(im.select("SCL").remap(MAU, [0] * len(MAU), 1)))
            .median())
    vis = comp.visualize(bands=["B4", "B3", "B2"], min=0, max=3000)
    agua = comp.normalizedDifference(["B3", "B11"]).gt(0).selfMask()
    vis = vis.blend(agua.visualize(palette=["1f77b4"], opacity=0.8))
    if outros:
        vis = vis.blend(ee.Image().paint(ee.FeatureCollection(outros), 0, 2).visualize(palette=["FF3030"]))
    vis = vis.blend(ee.Image().paint(ee.FeatureCollection([ee.Feature(env)]), 0, 2).visualize(palette=["FFFF00"]))
    vis = vis.blend(ee.Image().paint(ee.FeatureCollection([ee.Feature(ponto.buffer(120))]), 1)
                    .visualize(palette=["FF00FF"]))
    url = vis.getThumbURL({"region": caixa, "dimensions": 1000, "format": "png"})
    urllib.request.urlretrieve(url, caminho)


resumo = []
for alb in CFG["albufeiras"]:
    if args and alb["id"] not in args:
        continue
    pasta = os.path.join("dados", alb["id"])
    os.makedirs(pasta, exist_ok=True)
    destino = os.path.join(pasta, "envolvente.geojson")
    ponto = ee.Geometry.Point([alb["lon"], alb["lat"]])
    print("%s" % alb["nome"], flush=True)

    outros = []
    if os.path.exists(destino) and not refazer:
        print("  já existe (usar --refazer para recalcular)")
    elif alb["id"] == "alqueva" and os.path.exists("alqueva_envolvente.geojson"):
        shutil.copy("alqueva_envolvente.geojson", destino)       # a do passo 11
        print("  copiada de alqueva_envolvente.geojson (passo 11)")
    else:
        if "caixa" in alb:
            caixa = ee.Geometry.Rectangle(alb["caixa"])
        else:
            dx, dy = MEIA_CAIXA
            caixa = ee.Geometry.Rectangle([alb["lon"] - dx, alb["lat"] - dy, alb["lon"] + dx, alb["lat"] + dy])
        jrc = JRC.selfMask()
        for ex in alb.get("excluir", []):                 # zonas a retirar (ex.: rio a jusante de um açude)
            jrc = jrc.updateMask(ee.Image(1).paint(ee.FeatureCollection(
                [ee.Feature(ee.Geometry.Rectangle(ex))]), 0))
        vec = jrc.reduceToVectors(
            geometry=caixa, scale=30, geometryType="polygon", eightConnected=True,
            labelProperty="agua", maxPixels=1e10)
        vec = vec.filterBounds(ponto.buffer(RAIO_ESTACAO))
        vec = vec.map(lambda f: f.set("area_ha", f.geometry().area(30).divide(1e4)))
        cands = vec.sort("area_ha", False).limit(5)
        info = cands.getInfo()["features"]
        if not info:
            print("  !! nenhum polígono de água do JRC a menos de %d m da estação — rever lat/lon" % RAIO_ESTACAO)
            continue
        for i, f in enumerate(info):
            print("  candidato %d: %8.1f ha%s" % (i + 1, f["properties"]["area_ha"], "  <- escolhido" if i == 0 else ""))
        escolhido = ee.Feature(cands.first()).geometry()
        outros = [ee.Feature(ee.Geometry(f["geometry"])) for f in info[1:]]
        geo = escolhido.buffer(30, 10).simplify(30).getInfo()
        guardar_geojson(destino, geo, "Envolvente de %s (JRC max_extent + 30 m)" % alb["nome"])

    with open(destino, encoding="utf-8") as fh:
        env = ee.Geometry(json.load(fh)["features"][0]["geometry"], None, False)
    area_ha = env.area(10, PROJ).divide(1e4).getInfo()
    perim_km = env.perimeter(10, PROJ).divide(1e3).getInfo()
    print("  envolvente: %.0f ha | perímetro %.1f km | margem por área: %.2f km/km²"
          % (area_ha, perim_km, perim_km / (area_ha / 100)), flush=True)
    resumo.append([alb["id"], alb["nome"], round(area_ha, 1), round(perim_km, 2),
                   round(perim_km / (area_ha / 100), 3)])
    try:
        imagem(alb, env, outros, ponto, os.path.join(pasta, "envolvente.png"))
        print("  imagem: %s" % os.path.join(pasta, "envolvente.png"))
    except Exception as e:
        print("  (imagem falhou: %s)" % str(e)[:80])

if resumo:
    caminho = os.path.join("dados", "resumo_envolventes.csv")
    antigos = {}
    if os.path.exists(caminho):
        with open(caminho, encoding="utf-8") as fh:
            antigos = {r[0]: r for r in list(csv.reader(fh))[1:]}
    for r in resumo:
        antigos[r[0]] = r
    with open(caminho, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "nome", "area_envolvente_ha", "perimetro_km", "margem_km_por_km2"])
        w.writerows(antigos.values())
    print("\nGuardado: dados/<albufeira>/envolvente.geojson + envolvente.png, dados/resumo_envolventes.csv")
