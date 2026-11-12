"""
Monsaraz em Mudança — passo 20 (Fase 3)
O MÉTODO NOVO (NDWI a 10 m + filtro do núcleo) CONTRA O ORIGINAL (MNDWI a 20 m), EM TODAS AS ALBUFEIRAS.

Para cada albufeira e cada método:
  0. Validação da cota oficial (não usa o satélite): os sensores automáticos de cota têm falhas
     (derivas de semanas, saltos de vários metros). Uma data só entra se a cota for coerente com
     o volume oficial dos mesmos dias (cota vs V⁻¹(volume)) ou, sem volume por perto, com as
     cotas dos dias vizinhos (±2 dias, tolerância 1 m).
  1. Filtro automático (não usa a cota): fica a data se o núcleo de água aparecer como água em
     ≥ 98 %. Núcleo do passo 17b (Sentinel-2 no nível mais baixo, ≥ 2 ha) quando está completo;
     senão o do passo 17 (JRC, ≥ 5 ha). Sem núcleo, não há filtro.
     Uso: python 20_analise_ndwi10.py            (núcleo 17b)
          python 20_analise_ndwi10.py --nucleo-jrc (núcleo 17, para comparar)
  2. Curva área–cota A(h) de 2.º grau: r, dispersão.
  3. Área do satélite vs área implícita no volume oficial (dV/dh), na gama comum de cotas,
     e se a curva oficial é coerente (a área implícita cresce com a cota).
  4. Volume a partir de UMA imagem (área -> cota pela curva invertida -> volume pela V(h) oficial),
     validado deixando um ano de fora, com a curva ajustada SEM rejeição pela cota
     (só o filtro do núcleo): é o cenário de um observatório automático.

Não usa o Earth Engine. Requer dados/<id>/area_cota_ndwi10.csv (passo 19), qualidade.csv (17),
snirh_diario.csv (13) e dados/resumo_envolventes.csv (14).
Produz dados/resumo_fase3_ndwi10.csv e .json, e dados/<id>/cota_validada.csv
"""

import csv
import datetime as dt
import json
import os
from statistics import mean, median, pstdev

K_MAD = 3.5
NUCLEO_MIN_HA = 5          # núcleo do JRC (passo 17)
NUCLEO_MIN_HA_S2 = 2       # núcleo do Sentinel-2 (passo 17b)
import sys
NUCLEO = "jrc" if "--nucleo-jrc" in sys.argv else "s2"
NUCLEO_LIMIAR = 98


def polyfit(xs, ys, grau):
    n = grau + 1
    M = [[sum(x ** (i + j) for x in xs) for j in range(n)] + [sum(y * x ** i for x, y in zip(xs, ys))]
         for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def polyval(c, x):
    return sum(ci * x ** i for i, ci in enumerate(c))


def polyder(c):
    return [i * ci for i, ci in enumerate(c)][1:]


def ajuste_robusto(xs, ys, grau, x0):
    ok = [True] * len(xs)
    for _ in range(20):
        c = polyfit([x - x0 for x, k in zip(xs, ok) if k], [y for y, k in zip(ys, ok) if k], grau)
        res = [y - polyval(c, x - x0) for x, y in zip(xs, ys)]
        r_ok = [r for r, k in zip(res, ok) if k]
        med = median(r_ok)
        mad = 1.4826 * median(abs(r - med) for r in r_ok)
        novo = [abs(r - med) <= K_MAD * mad for r in res]
        if novo == ok:
            break
        ok = novo
    return c, res, ok


def pearson(xs, ys):
    mx, my = mean(xs), mean(ys)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    return sxy / (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5


def quantil(v, p):
    v = sorted(v)
    return v[int(p * (len(v) - 1))]


def rmse(e):
    return (sum(x * x for x in e) / len(e)) ** 0.5


def inverter(c, h0, a, lo, hi):
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if polyval(c, m - h0) < a else (lo, m)
    return (lo + hi) / 2



def validar_cotas(pasta, tem_volume):
    """Devolve {data: (estado, motivo)} para cada dia com cota. estado: True/False/None (não verificável)."""
    with open(os.path.join(pasta, "snirh_diario.csv"), encoding="utf-8") as fh:
        linhas = [r for r in csv.DictReader(fh) if r["cota_m"]]
    cota = {dt.date.fromisoformat(r["data"]): float(r["cota_m"]) for r in linhas}
    h_vol, tol = {}, None
    if tem_volume:
        pares = [(float(r["cota_m"]), float(r["volume_dam3"]) / 1000, dt.date.fromisoformat(r["data"]))
                 for r in linhas if r["volume_dam3"]]
        if len(pares) >= 15:
            Ho, Vo = [p[0] for p in pares], [p[1] for p in pares]
            hv0 = round(median(Ho), 2)
            cV, _, okV = ajuste_robusto(Ho, Vo, 3, hv0)
            Hok = [h for h, k in zip(Ho, okV) if k]
            lo, hi = min(Hok) - 2, max(Hok) + 2

            def inv(v):
                a, b = lo, hi
                for _ in range(50):
                    m = (a + b) / 2
                    a, b = (m, b) if polyval(cV, m - hv0) < v else (a, m)
                return (a + b) / 2
            h_vol = {d: inv(v) for _, v, d in pares}
            dif = [h - h_vol[d] for h, _, d in pares]
            med = median(dif)
            tol = max(0.30, K_MAD * 1.4826 * median(abs(x - med) for x in dif))
    estado = {}
    for d, c in cota.items():
        perto = [h_vol[d + dt.timedelta(days=k)] for k in (0, -1, 1, -2, 2, -3, 3) if d + dt.timedelta(days=k) in h_vol]
        if perto:
            estado[d] = (abs(c - perto[0]) <= tol, "volume")
            continue
        viz = [cota[d + dt.timedelta(days=k)] for k in (-2, -1, 1, 2) if d + dt.timedelta(days=k) in cota]
        if len(viz) >= 2:
            estado[d] = (abs(c - median(viz)) <= 1.0, "vizinhos")
        else:
            estado[d] = (None, "sem verificação")
    return estado, tol


with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)
with open(os.path.join("dados", "resumo_envolventes.csv"), encoding="utf-8") as fh:
    ENV = {r["id"]: r for r in csv.DictReader(fh)}


def analisar(alb, coluna):
    pasta = os.path.join("dados", alb["id"])
    with open(os.path.join(pasta, "area_cota_ndwi10.csv"), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    # núcleo: o do passo 17b (Sentinel-2, nível mínimo) se estiver completo; senão o do passo 17 (JRC)
    fonte, nucleo_pct = "jrc", {r["data"]: r["nucleo_agua_pct"] for r in rows}
    with open(os.path.join(pasta, "qualidade.csv"), encoding="utf-8") as fh:
        q = list(csv.DictReader(fh))
    nucleo_ha = float(q[0]["nucleo_ha"]) if q else 0
    minimo = NUCLEO_MIN_HA
    f_s2 = os.path.join(pasta, "qualidade_s2.csv")
    if NUCLEO == "s2" and os.path.exists(f_s2):
        with open(f_s2, encoding="utf-8") as fh:
            q2 = {r["data"]: r for r in csv.DictReader(fh)}
        if all(r["data"] in q2 for r in rows):
            fonte, minimo = "s2", NUCLEO_MIN_HA_S2
            nucleo_ha = float(next(iter(q2.values()))["nucleo_ha"])
            nucleo_pct = {d: r["nucleo_agua_pct"] for d, r in q2.items()}
        else:
            fonte = "jrc (17b incompleto)"
    elif NUCLEO == "s2":
        fonte = "jrc (sem núcleo S2)"     # 17b não correu, ou a albufeira é pequena demais no mínimo
    for r in rows:
        r["nucleo_agua_pct"] = nucleo_pct.get(r["data"], "")
    filtra = nucleo_ha >= minimo
    estado, tol = VALID[alb["id"]]

    def cota_ok(r):          # a cota usada veio do dia data + dif_dias
        d = dt.date.fromisoformat(r["data"]) + dt.timedelta(days=int(r["dif_dias"]))
        return estado.get(d, (None, ""))[0] is not False
    n_total = len(rows)
    rows = [r for r in rows if cota_ok(r)]
    usar = [r for r in rows if not filtra or (r["nucleo_agua_pct"] and float(r["nucleo_agua_pct"]) >= NUCLEO_LIMIAR)]
    D = [r["data"] for r in usar]
    H = [float(r["cota_m"]) for r in usar]
    A = [float(r[coluna]) for r in usar]
    h0 = round(mean(H), 2)
    out = {"fonte_nucleo": fonte, "datas": n_total, "datas_cota_rejeitada": n_total - len(rows), "datas_apos_filtro": len(usar),
           "filtro_nucleo": filtra, "nucleo_ha": round(nucleo_ha, 1),
           "tolerancia_cota_m": None if tol is None else round(tol, 2),
           "cota_min_todas_m": min(float(r["cota_m"]) for r in rows), "cota_min_apos_filtro_m": min(H)}

    cA, res, ok = ajuste_robusto(H, A, 2, h0)
    Hc = [h for h, k in zip(H, ok) if k]
    Ac = [a for a, k in zip(A, ok) if k]
    disp = pstdev([r for r, k in zip(res, ok) if k])
    out.update({"r_apos_filtro": round(pearson(H, A), 3), "anomalas_restantes": len(H) - len(Hc),
                "r_sem_anomalas": round(pearson(Hc, Ac), 3),
                "dispersao_pct": round(100 * disp / mean(Ac), 1)})

    if not alb.get("volume"):
        return out
    with open(os.path.join(pasta, "snirh_diario.csv"), encoding="utf-8") as fh:
        of = [r for r in csv.DictReader(fh) if r["cota_m"] and r["volume_dam3"]
              and estado.get(dt.date.fromisoformat(r["data"]), (None, ""))[0] is not False]
    Ho = [float(r["cota_m"]) for r in of]
    Vo = [float(r["volume_dam3"]) / 1000 for r in of]
    hv0 = round(median(Ho), 2)
    cV, _, okV = ajuste_robusto(Ho, Vo, 3, hv0)
    V = lambda h: polyval(cV, h - hv0)
    dV = polyder(cV)
    Hok = [h for h, k in zip(Ho, okV) if k]

    lo = max(quantil(Hc, 0.10), quantil(Hok, 0.02))
    hi = min(quantil(Hc, 0.90), quantil(Hok, 0.98))
    comp = []
    for i in range(8):
        h = lo + (hi - lo) * i / 7
        comp.append((polyval(cA, h - h0), polyval(dV, h - hv0) * 100))
    por_data = [a / (polyval(dV, h - hv0) * 100) for h, a, k in zip(H, A, ok) if k and lo <= h <= hi]
    out["razao_area_satelite_oficial"] = round(median(por_data), 3) if por_data else None
    out["curva_oficial_consistente"] = all(b2 > b1 for (_, b1), (_, b2) in zip(comp, comp[1:]))

    # volume a partir de uma imagem, um ano de fora, curva SEM rejeição pela cota
    hmin, hmax = min(Hok) - 1, max(Hok) + 1
    erros, erros_mad = [], []
    for ano in sorted(set(d[:4] for d in D)):
        tr = [i for i in range(len(H)) if D[i][:4] != ano]
        te = [i for i in range(len(H)) if D[i][:4] == ano]
        if len(tr) < 10:
            continue
        c = polyfit([H[i] - h0 for i in tr], [A[i] for i in tr], 2)
        erros += [V(inverter(c, h0, A[i], hmin, hmax)) - V(H[i]) for i in te]
        cr, rr, okr = ajuste_robusto([H[i] for i in tr], [A[i] for i in tr], 2, h0)
        r_ok = [x for x, k in zip(rr, okr) if k]
        med, mad = median(r_ok), 1.4826 * median(abs(x - median(r_ok)) for x in r_ok)
        erros_mad += [V(inverter(cr, h0, A[i], hmin, hmax)) - V(H[i]) for i in te
                      if abs(A[i] - polyval(cr, H[i] - h0) - med) <= K_MAD * mad]
    vm = mean(V(h) for h in H)
    amp = V(quantil(H, 0.98)) - V(quantil(H, 0.02))
    out.update({"vol_rmse_hm3": round(rmse(erros), 2), "vol_rmse_pct_medio": round(100 * rmse(erros) / vm, 1),
                "vol_rmse_pct_variacao": round(100 * rmse(erros) / amp, 1) if amp > 0 else None,
                "vol_vies_hm3": round(mean(erros), 2), "vol_erro_mediano_pct_medio": round(100 * median(abs(e) for e in erros) / vm, 1),
                "vol_rmse_pct_medio_com_mad": round(100 * rmse(erros_mad) / vm, 1), "volume_medio_hm3": round(vm, 1)})
    return out


VALID = {}
resultados = []
for alb in CFG["albufeiras"]:
    if not os.path.exists(os.path.join("dados", alb["id"], "area_cota_ndwi10.csv")):
        continue
    VALID[alb["id"]] = validar_cotas(os.path.join("dados", alb["id"]), alb.get("volume"))
    with open(os.path.join("dados", alb["id"], "cota_validada.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["data", "cota_coerente", "verificada_por"])
        for d, (e, motivo) in sorted(VALID[alb["id"]][0].items()):
            w.writerow([d.isoformat(), "" if e is None else ("sim" if e else "nao"), motivo])
    r = {"id": alb["id"], "nome": alb["nome"], "envolvente_ha": float(ENV[alb["id"]]["area_envolvente_ha"]),
         "margem_km_por_km2": float(ENV[alb["id"]]["margem_km_por_km2"]),
         "mndwi20": analisar(alb, "agua_mndwi20_ha"), "ndwi10": analisar(alb, "agua_ha")}
    resultados.append(r)

print("Erro do volume: um ano de fora, só com o filtro do núcleo; em % do volume médio e (entre parênteses) em % da variação de volume do período.")
print("%-11s %5s %6s | %-19s | %-15s | %-31s" % ("", "", "cota", "satélite / oficial", " r (só filtro)", "       erro do volume"))
print("%-11s %5s %6s | %8s %9s | %7s %7s | %15s %15s" % ("albufeira", "datas", "rejeit", "MNDWI20", "NDWI10", "MNDWI20", "NDWI10", "MNDWI20", "NDWI10"))
for r in resultados:
    m, n = r["mndwi20"], r["ndwi10"]
    f = lambda d, k, fmt: (fmt % d[k]) if d.get(k) is not None else "—"
    ev = lambda d: ("%.1f%% (%s%%)" % (d["vol_rmse_pct_medio"], f(d, "vol_rmse_pct_variacao", "%.0f"))) if d.get("vol_rmse_pct_medio") is not None else "—"
    notas = [] if n["fonte_nucleo"] in ("s2", "jrc") else ["núcleo: " + n["fonte_nucleo"]]
    if not n.get("curva_oficial_consistente", True):
        notas.append("curva oficial incoerente")
    if not n["filtro_nucleo"]:
        notas.append("sem filtro (sem núcleo)")
    elif n["cota_min_apos_filtro_m"] - n["cota_min_todas_m"] > 0.5:
        notas.append("filtro corta as cotas mais baixas (%.1f -> %.1f m)" % (n["cota_min_todas_m"], n["cota_min_apos_filtro_m"]))
    print("%-11s %5d %6d | %8s %9s | %7s %7s | %15s %15s  %s" % (
        r["nome"], n["datas_apos_filtro"], n["datas_cota_rejeitada"], f(m, "razao_area_satelite_oficial", "%.2f"),
        f(n, "razao_area_satelite_oficial", "%.2f"), f(m, "r_apos_filtro", "%.2f"), f(n, "r_apos_filtro", "%.2f"),
        ev(m), ev(n), "; ".join(notas)))

campos = ["id", "nome", "envolvente_ha", "margem_km_por_km2"]
with open(os.path.join("dados", "resumo_fase3_ndwi10.csv"), "w", newline="", encoding="utf-8") as fh:
    chaves = sorted({k for r in resultados for k in r["ndwi10"]})
    w = csv.writer(fh)
    w.writerow(campos + ["mndwi20_" + k for k in chaves] + ["ndwi10_" + k for k in chaves])
    for r in resultados:
        w.writerow([r[c] for c in campos] + [r["mndwi20"].get(k, "") for k in chaves] + [r["ndwi10"].get(k, "") for k in chaves])
with open(os.path.join("dados", "resumo_fase3_ndwi10.json"), "w", encoding="utf-8") as fh:
    json.dump({"albufeiras": resultados,
               "metodo": "NDWI (B3,B8) > -0,1 a 10 m, limiar calibrado só no Alqueva e no Caia; filtro: núcleo de água permanente (JRC) visto como água em >= 98 %.",
               "fontes": "Área: Sentinel-2 (análise própria). Cota e volume: SNIRH (APA, I.P.). Envolventes e núcleo: JRC Global Surface Water."},
              fh, ensure_ascii=False, indent=2)
print("\nGuardado: dados/resumo_fase3_ndwi10.csv e .json")
