"""
Monsaraz em Mudança — passo 05
Comparacao temporal entre anos, no MESMO trimestre (maçãs com maçãs).

Produz, para cada ano escolhido:
  - agua_AAAA_T3.geojson  -> poligonos da agua delimitada (dados da app)
Produz, no conjunto:
  - comparacao_T3.json    -> areas e diferencas (indicadores da app)
  - mudanca_2020_2025.png -> mapa de mudança (estavel/ganho/perda)
"""

import ee
import json
import urllib.request

ee.Initialize(project="projetobeta-513108")

with open("monsaraz.geojson", encoding="utf-8") as fh:
    geom = json.load(fh)["features"][0]["geometry"]
aoi = ee.Geometry(geom, None, False)

PROJ = "EPSG:3763"
ANOS = list(range(2017, 2026))            # 9 anos: 2017..2025
Q = 3                                     # trimestre a comparar (1..4)
TRIMESTRES = {1: ("01-01", "03-31"), 2: ("04-01", "06-30"),
              3: ("07-01", "09-30"), 4: ("10-01", "12-31")}


def mascara_nuvens(img):
    scl = img.select("SCL")
    mau = (scl.eq(0).Or(scl.eq(1)).Or(scl.eq(3))
           .Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)))
    return img.updateMask(mau.Not()).select(["B2", "B3", "B4", "B11"])


def composicao(ano, q):
    ini, fim = TRIMESTRES[q]
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
           .filterBounds(aoi).filterDate(f"{ano}-{ini}", f"{ano}-{fim}")
           .map(mascara_nuvens))
    return col.median().clip(aoi)


def agua_mask(comp):
    return comp.normalizedDifference(["B3", "B11"]).gt(0).unmask(0).clip(aoi)


def area_ha(mask):
    m2 = (mask.selfMask().multiply(ee.Image.pixelArea())
          .reduceRegion(ee.Reducer.sum(), aoi, 10, PROJ, maxPixels=1e9)
          .values().get(0))
    return ee.Number(m2).divide(10000)


# --- 1. Por ano: area + vetorizacao -> geojson -----------------------
areas = {}
masks = {}
comps = {}
for ano in ANOS:
    comp = composicao(ano, Q)
    m = agua_mask(comp)
    comps[ano] = comp
    masks[ano] = m
    a = area_ha(m).getInfo()
    areas[ano] = round(a, 1)
    print("%d-T%d: %.1f ha" % (ano, Q, a))

    # vetorizar: dissolver + simplificar (tolerancia 10 m)
    vec = m.selfMask().reduceToVectors(
        geometry=aoi, scale=10, geometryType="polygon",
        eightConnected=True, crs="EPSG:4326", maxPixels=1e9, bestEffort=True)
    g = vec.geometry().simplify(10).getInfo()
    feat = {"type": "Feature",
            "properties": {"ano": ano, "trimestre": "T%d" % Q,
                           "area_ha": round(a, 1), "fonte": "Sentinel-2 / MNDWI"},
            "geometry": g}
    fc = {"type": "FeatureCollection", "features": [feat]}
    nome = "agua_%d_T%d.geojson" % (ano, Q)
    with open(nome, "w", encoding="utf-8") as fh:
        json.dump(fc, fh, ensure_ascii=False)
    print("   -> %s" % nome)

# --- 2. Indicadores: diferencas --------------------------------------
def dif(a, b):
    d = areas[b] - areas[a]
    return {"de": a, "para": b, "delta_ha": round(d, 1),
            "delta_pct": round(100 * d / areas[a], 1)}

comp_json = {"trimestre": "T%d" % Q, "areas_ha": areas,
             "diferencas": [dif(ANOS[0], ANOS[-1]),
                            dif(ANOS[0], ANOS[1]), dif(ANOS[1], ANOS[-1])]}
with open("comparacao_T%d.json" % Q, "w", encoding="utf-8") as fh:
    json.dump(comp_json, fh, ensure_ascii=False, indent=2)

print("\nDiferencas:")
for d in comp_json["diferencas"]:
    print("  %d -> %d: %+.1f ha (%+.1f%%)" % (d["de"], d["para"], d["delta_ha"], d["delta_pct"]))

# --- 3. Mapa de mudança (2020 seco vs 2025 cheio) --------------------
a0, a1 = 2020, 2025
w0, w1 = masks[a0], masks[a1]
estavel = w0.And(w1)
ganho = w1.And(w0.Not())
perda = w0.And(w1.Not())

rgb = comps[a1].visualize(bands=["B4", "B3", "B2"], min=0, max=3000)
vis = (rgb
       .blend(estavel.selfMask().visualize(palette=["1f77b4"]))   # azul
       .blend(ganho.selfMask().visualize(palette=["2ca02c"]))     # verde
       .blend(perda.selfMask().visualize(palette=["d62728"])))    # vermelho
contorno = ee.Image().paint(ee.FeatureCollection([ee.Feature(aoi)]), 0, 2)
vis = vis.blend(contorno.visualize(palette=["FFFF00"]))

url = vis.getThumbURL({"region": aoi.bounds(), "dimensions": 1024, "format": "png"})
urllib.request.urlretrieve(url, "mudanca_%d_%d.png" % (a0, a1))
print("\nGuardado: mudanca_%d_%d.png (azul=estavel, verde=ganho, vermelho=perda)" % (a0, a1))
