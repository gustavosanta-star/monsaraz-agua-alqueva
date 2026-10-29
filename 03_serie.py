"""
Monsaraz em Mudança — passo 03
Serie temporal da area de agua: a mesma logica, agora PARAMETRIZADA por ano.

Objetivo de engenharia: acrescentar uma data = mais uma volta do ciclo,
sem repetir trabalho manual. Grava os resultados em area_serie.csv.
"""

import ee
import json
import csv

ee.Initialize(project="projetobeta-513108")

# Limite da freguesia (CAOP)
with open("monsaraz.geojson", encoding="utf-8") as fh:
    geom = json.load(fh)["features"][0]["geometry"]
aoi = ee.Geometry(geom, None, False)


def agua_no_verao(ano):
    """Devolve (data, nuvens_%, area_ha) para a imagem de verao menos nublada do ano."""
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
           .filterBounds(aoi)
           .filterDate(f"{ano}-06-01", f"{ano}-09-30")
           .sort("CLOUDY_PIXEL_PERCENTAGE"))

    img = ee.Image(col.first()).clip(aoi)

    mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI")
    agua = mndwi.gt(0)

    area_m2 = (agua.multiply(ee.Image.pixelArea())
               .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                             scale=10, crs="EPSG:3763", maxPixels=1e9)
               .get("MNDWI"))

    info = ee.Dictionary({
        "data": ee.Date(img.get("system:time_start")).format("YYYY-MM-dd"),
        "nuvens": img.get("CLOUDY_PIXEL_PERCENTAGE"),
        "area_ha": ee.Number(area_m2).divide(10000),
    }).getInfo()
    return info


linhas = []
print("ano | data        | nuvens | area (ha)")
print("----+-------------+--------+----------")
for ano in range(2017, 2026):
    try:
        r = agua_no_verao(ano)
        print("%4d| %s |  %4.1f%% | %8.1f" % (ano, r["data"], r["nuvens"], r["area_ha"]))
        linhas.append([ano, r["data"], round(r["nuvens"], 1), round(r["area_ha"], 1)])
    except Exception as e:
        print("%4d| sem imagem adequada (%s)" % (ano, e))

with open("area_serie.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["ano", "data", "nuvens_pct", "area_ha"])
    w.writerows(linhas)

print("\nGuardado: area_serie.csv (%d anos)" % len(linhas))
