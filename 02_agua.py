"""
Monsaraz em Mudança — passo 02
Limite oficial da freguesia (CAOP) + deteção de água com MNDWI.

O que este script faz:
  1. Lê o limite da freguesia a partir de monsaraz.geojson (tem de estar na mesma pasta).
  2. Escolhe a imagem Sentinel-2 menos nublada de um verão.
  3. Calcula o índice de água MNDWI e separa água (MNDWI > 0).
  4. Mede a área de água dentro da freguesia, em hectares.
  5. Grava um PNG: imagem a cores + água a azul + contorno da freguesia a amarelo.
"""

import ee
import json
import urllib.request

ee.Initialize(project="projetobeta-513108")

# --- 1. Limite da freguesia (lido do ficheiro local) -----------------
with open("monsaraz.geojson", encoding="utf-8") as fh:
    gj = json.load(fh)
geom = gj["features"][0]["geometry"]
# geodesic=False: arestas planas, correto para um limite administrativo
aoi = ee.Geometry(geom, None, False)

# --- 2. Melhor imagem Sentinel-2 do verao ----------------------------
col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
       .filterBounds(aoi)
       .filterDate("2024-06-01", "2024-09-30")
       .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 10))
       .sort("CLOUDY_PIXEL_PERCENTAGE"))

img = col.first().clip(aoi)
data = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd").getInfo()
print("Imagem escolhida, data:", data)

# --- 3. MNDWI = (Verde - SWIR1) / (Verde + SWIR1) --------------------
# Sentinel-2: B3 = Verde, B11 = SWIR1. Agua tem MNDWI alto.
mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI")
agua = mndwi.gt(0)                      # limiar inicial: 0

# --- 4. Area de agua dentro da freguesia, em hectares ----------------
area_m2 = (agua.multiply(ee.Image.pixelArea())
           .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                         scale=10, crs="EPSG:3763", maxPixels=1e9)
           .get("MNDWI"))
area_ha = ee.Number(area_m2).divide(10000).getInfo()
print("Area de agua estimada: %.1f ha" % area_ha)

# --- 5. Imagem final: RGB + agua a azul + contorno a amarelo ---------
rgb = img.visualize(bands=["B4", "B3", "B2"], min=0, max=3000)
agua_vis = agua.selfMask().visualize(palette=["00BFFF"])        # azul
base = rgb.blend(agua_vis)
contorno = ee.Image().paint(ee.FeatureCollection([ee.Feature(aoi)]), 0, 2)
final = base.blend(contorno.visualize(palette=["FFFF00"]))      # amarelo

url = final.getThumbURL({"region": aoi.bounds(), "dimensions": 1024, "format": "png"})
urllib.request.urlretrieve(url, "monsaraz_agua.png")
print("Guardado: monsaraz_agua.png")
