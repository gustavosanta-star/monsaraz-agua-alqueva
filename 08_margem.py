"""
Monsaraz em Mudança — passo 08
ERRO DE MARGEM: termo de incerteza baseado no comprimento da margem.

Porquê: o MNDWI usa a banda B11, que é nativa a 20 m. Ao longo de uma margem
muito recortada, os pixéis mistos (meio água, meio terra) podem pesar tanto ou
mais do que a escolha do limiar (o ± do passo 07).

Como:
  1. Para cada ano (T3), deteta a água numa zona um pouco maior do que a
     freguesia (buffer de 500 m), para que o limite administrativo não crie
     margens artificiais no meio do lago.
  2. Passa a máscara para a grelha de 20 m (resolução nativa da B11).
  3. Conta os pixéis de margem dentro da freguesia: pixéis de água com pelo
     menos um vizinho (4-ligação) que não é água.
  4. Comprimento da margem L ≈ n.º de pixéis de margem × 20 m.
     (Em margens diagonais, esta contagem em escada sobrestima L; os limites
     abaixo ficam por isso do lado conservador.)

Dois números por ano:
  - limite_sup_ha: enviesamento SISTEMÁTICO de meio pixel (10 m) ao longo de
    toda a margem = L × 10 m. É o pior caso: todos os pixéis mistos errados no
    mesmo sentido.
  - aleatorio_ha: erros INDEPENDENTES que se anulam parcialmente. Cada pixel de
    margem tem uma fração de água desconhecida (desvio-padrão ≈ 0,29 de um
    pixel de 400 m²); o total cresce com a raiz do n.º de pixéis.
O erro real deverá estar entre os dois. O ± do limiar (passo 07) é mostrado
ao lado para comparação.
"""

import ee
import json
import csv
import os

ee.Initialize(project="projetobeta-513108")

with open("monsaraz.geojson", encoding="utf-8") as fh:
    geom = json.load(fh)["features"][0]["geometry"]
aoi = ee.Geometry(geom, None, False)
zona = aoi.buffer(500)

GRELHA = ee.Projection("EPSG:3763").atScale(20)      # 20 m, como a B11
S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
PIX_M = 20.0
DESVIO_FRACAO = 0.29                                  # desvio-padrão de uma fração uniforme 0–1


def mascara(img):
    scl = img.select("SCL")
    mau = (scl.eq(0).Or(scl.eq(1)).Or(scl.eq(3)).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)))
    return img.updateMask(mau.Not()).select(["B3", "B11"])


def pixeis_margem(ano):
    comp = (S2.filterBounds(zona).filterDate(f"{ano}-07-01", f"{ano}-09-30")
            .map(mascara).median())
    agua = (comp.normalizedDifference(["B3", "B11"]).gt(0)
            .unmask(0).clip(zona).reproject(GRELHA))
    vizinhos_todos_agua = agua.focal_min(kernel=ee.Kernel.plus(1, "pixels")).reproject(GRELHA)
    margem = agua.And(vizinhos_todos_agua.Not()).rename("margem")
    n = (margem.selfMask()
         .reduceRegion(ee.Reducer.count(), aoi, scale=PIX_M, crs="EPSG:3763", maxPixels=1e9)
         .get("margem"))
    return ee.Number(n).getInfo()


# ± do limiar (passo 07), se existir, só para comparação
limiar = {}
if os.path.exists("incerteza.csv"):
    with open("incerteza.csv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            limiar[int(r["ano"])] = float(r["incerteza_ha"])

linhas = []
print("ano | margem (km) | sistemático ± (ha) | aleatório ± (ha) | limiar ± (ha)")
print("----+-------------+--------------------+------------------+--------------")
for ano in range(2017, 2026):
    n = pixeis_margem(ano)
    L_m = n * PIX_M
    sup_ha = L_m * (PIX_M / 2) / 10000
    ale_ha = (n ** 0.5) * DESVIO_FRACAO * PIX_M * PIX_M / 10000
    lim = limiar.get(ano)
    print("%4d| %11.1f | %18.1f | %16.1f | %12s" %
          (ano, L_m / 1000, sup_ha, ale_ha, ("%.1f" % lim) if lim is not None else "—"))
    linhas.append([ano, round(L_m / 1000, 1), round(sup_ha, 1), round(ale_ha, 1)])

with open("margem.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["ano", "margem_km", "sistematico_ha", "aleatorio_ha"])
    w.writerows(linhas)

media_sup = sum(l[2] for l in linhas) / len(linhas)
media_km = sum(l[1] for l in linhas) / len(linhas)
print("\nMargem média: %.1f km dentro da freguesia" % media_km)
print("Limite superior médio (meio pixel de 20 m, sistemático): ±%.1f ha" % media_sup)
print("Guardado: margem.csv")
