"""
Monsaraz em Mudança — passo 16 (Fase 3)
CURVAS ÁREA–COTA E VOLUME PARA TODAS AS ALBUFEIRAS — e o teste da hipótese da margem.

É o passo 12 generalizado (mesmo método), aplicado a cada albufeira do albufeiras.json:
  Fase 1  A(h) de 2.º grau com rejeição robusta (3,5 × MAD) das datas anómalas;
  Fase 2  V(h) oficial de 3.º grau (mesma rejeição); área do satélite vs área implícita
          no volume oficial (dV/dh); ΔV = ∫A dh; série de volume a partir de UM volume
          de referência.
Datas utilizáveis: lacuna < 0,2 % da envolvente (o mesmo critério do Alqueva: 50 ha em 24 596).

Hipótese em teste: o satélite subestima a área por causa dos pixels mistos da margem,
por isso a subestimativa deve crescer com a margem por unidade de área (km de margem por km²).

Não usa o Earth Engine. Requer dados/<id>/area_cota.csv (passo 15), snirh_diario.csv
(passo 13) e dados/resumo_envolventes.csv (passo 14).

Produz:
  - dados/<id>/area_cota_ajuste.csv, dados/<id>/volume_comparacao.csv, dados/<id>/resumo.json
  - dados/resumo_fase3.csv e dados/resumo_fase3.json (uma linha por albufeira + teste da hipótese)
"""

import csv
import json
import os
from statistics import mean, median, pstdev

LACUNA_FRAC = 50 / 24596
K_MAD = 3.5
MIN_DATAS = 15


# --- utilitários (sem numpy), iguais ao passo 12 ------------------------------
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


def polyint(c, a, b):
    return sum(ci / (i + 1) * (b ** (i + 1) - a ** (i + 1)) for i, ci in enumerate(c))


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
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / (sxx * syy) ** 0.5


def quantil(v, p):
    v = sorted(v)
    return v[int(p * (len(v) - 1))]


# --- por albufeira -----------------------------------------------------------
with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)
with open(os.path.join("dados", "resumo_envolventes.csv"), encoding="utf-8") as fh:
    ENV = {r["id"]: r for r in csv.DictReader(fh)}


def analisar(alb):
    pasta = os.path.join("dados", alb["id"])
    env_ha = float(ENV[alb["id"]]["area_envolvente_ha"])
    with open(os.path.join(pasta, "area_cota.csv"), encoding="utf-8") as fh:
        todas = list(csv.DictReader(fh))
    uteis = [r for r in todas if r["cota_m"] and float(r["lacuna_ha"]) < LACUNA_FRAC * env_ha]
    out = {"id": alb["id"], "nome": alb["nome"], "envolvente_ha": env_ha,
           "margem_km_por_km2": float(ENV[alb["id"]]["margem_km_por_km2"]),
           "datas_total": len(todas), "datas_utilizaveis": len(uteis)}
    if len(uteis) < MIN_DATAS:
        out["nota"] = "poucas datas utilizáveis (< %d)" % MIN_DATAS
        return out

    # Fase 1
    H = [float(r["cota_m"]) for r in uteis]
    A = [float(r["agua_ha"]) for r in uteis]
    h0 = round(mean(H), 2)
    cA, resA, okA = ajuste_robusto(H, A, 2, h0)
    Hc = [h for h, k in zip(H, okA) if k]
    Ac = [a for a, k in zip(A, okA) if k]
    res_c = [r for r, k in zip(resA, okA) if k]
    disp = pstdev(res_c)
    ss_tot = sum((a - mean(Ac)) ** 2 for a in Ac)
    out.update({"datas_na_curva": len(Hc), "datas_anomalas": len(H) - len(Hc),
                "pct_anomalas": round(100 * (len(H) - len(Hc)) / len(H), 1),
                "r_todas": round(pearson(H, A), 3), "r_sem_anomalas": round(pearson(Hc, Ac), 3),
                "r2_quadratico": round(1 - sum(r * r for r in res_c) / ss_tot, 3),
                "dispersao_ha": round(disp, 1), "dispersao_pct": round(100 * disp / mean(Ac), 1),
                "cota_min_m": min(Hc), "cota_max_m": max(Hc),
                "curva": {"h0": h0, "c": [round(x, 4) for x in cA]}})
    with open(os.path.join(pasta, "area_cota_ajuste.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["data", "cota_m", "agua_ha", "prevista_ha", "residuo_ha", "anomala"])
        for r, h, a, res, k in zip(uteis, H, A, resA, okA):
            w.writerow([r["data"], h, a, round(a - res, 2), round(res, 2), "nao" if k else "sim"])

    # Fase 2
    if not alb.get("volume"):
        out["nota"] = "sem volume no SNIRH: só a Fase 1"
        return out
    with open(os.path.join(pasta, "snirh_diario.csv"), encoding="utf-8") as fh:
        of = [r for r in csv.DictReader(fh) if r["cota_m"] and r["volume_dam3"]]
    Ho = [float(r["cota_m"]) for r in of]
    Vo = [float(r["volume_dam3"]) / 1000 for r in of]
    if len(Ho) < MIN_DATAS:
        out["nota"] = "poucos registos de volume"
        return out
    hv0 = round(median(Ho), 2)
    cV, resV, okV = ajuste_robusto(Ho, Vo, 3, hv0)
    dV = polyder(cV)
    Hok = [h for h, k in zip(Ho, okV) if k]

    # comparação de áreas só onde há satélite e volume oficial (percentis 10–90 da sobreposição)
    lo = max(quantil(Hc, 0.10), quantil(Hok, 0.02))
    hi = min(quantil(Hc, 0.90), quantil(Hok, 0.98))
    comp = []
    if hi > lo:
        for i in range(8):
            h = lo + (hi - lo) * i / 7
            a_sat = polyval(cA, h - h0)
            a_of = polyval(dV, h - hv0) * 100            # hm³/m -> ha
            if a_of > 0:
                comp.append({"cota_m": round(h, 2), "area_satelite_ha": round(a_sat, 1),
                             "area_oficial_ha": round(a_of, 1), "razao": round(a_sat / a_of, 3)})
    razao = mean(c["razao"] for c in comp) if comp else None
    # a área implícita no volume oficial tem de crescer com a cota; se não cresce, a curva
    # oficial V(h) é inconsistente nessa gama e a comparação de áreas não é fiável
    of_areas = [c["area_oficial_ha"] for c in comp]
    oficial_ok = all(b2 > a2 for a2, b2 in zip(of_areas, of_areas[1:]))
    # mesma comparação em volume: ΔV = ∫A dh (satélite) vs ΔV oficial, na gama comum
    razao_dv = None
    if hi > lo:
        dv_of = polyval(cV, hi - hv0) - polyval(cV, lo - hv0)
        if dv_of > 0:
            razao_dv = polyint(cA, lo - h0, hi - h0) / 100 / dv_of

    # série de volume a partir de um volume de referência (registo mais perto da mediana da cota)
    dentro = [i for i in range(len(Ho)) if okV[i] and min(Hc) <= Ho[i] <= max(Hc)]
    i_ref = min(dentro, key=lambda i: abs(Ho[i] - median(Hc)))
    v_sat = [Vo[i_ref] + polyint(cA, Ho[i_ref] - h0, h - h0) / 100 for h in Ho]
    erros = [vs - vo for i, (vs, vo) in enumerate(zip(v_sat, Vo)) if i in set(dentro)]
    rmse = (sum(e * e for e in erros) / len(erros)) ** 0.5
    v_dentro = [Vo[i] for i in dentro]
    amplitude_v = quantil(v_dentro, 0.98) - quantil(v_dentro, 0.02)
    out.update({"registos_volume": len(Ho), "registos_volume_anomalos": len(Ho) - sum(okV),
                "razao_area_satelite_oficial": None if razao is None else round(razao, 3),
                "razao_dv_satelite_oficial": None if razao_dv is None else round(razao_dv, 3),
                "curva_oficial_consistente": oficial_ok,
                "comparacao_por_cota": comp,
                "serie_registos": len(erros), "serie_referencia": of[i_ref]["data"],
                "serie_erro_medio_hm3": round(mean(erros), 2), "serie_rmse_hm3": round(rmse, 2),
                "serie_rmse_pct_volume_medio": round(100 * rmse / mean(v_dentro), 1),
                "serie_rmse_pct_variacao": round(100 * rmse / amplitude_v, 1) if amplitude_v else None})
    with open(os.path.join(pasta, "volume_comparacao.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["data", "cota_m", "volume_oficial_hm3", "volume_satelite_hm3", "anomalo_oficial", "na_gama_satelite"])
        for i, (r, h, vo, vs, k) in enumerate(zip(of, Ho, Vo, v_sat, okV)):
            w.writerow([r["data"], h, round(vo, 2), round(vs, 2), "nao" if k else "sim",
                        "sim" if i in set(dentro) else "nao"])
    return out


resultados = []
for alb in CFG["albufeiras"]:
    if not os.path.exists(os.path.join("dados", alb["id"], "area_cota.csv")):
        print("%-11s sem area_cota.csv (correr o passo 15)" % alb["nome"])
        continue
    r = analisar(alb)
    with open(os.path.join("dados", alb["id"], "resumo.json"), "w", encoding="utf-8") as fh:
        json.dump(r, fh, ensure_ascii=False, indent=2)
    resultados.append(r)

print("\n%-11s %8s %6s %5s %6s %6s %6s %7s %7s %7s" % (
    "albufeira", "área_ha", "margem", "datas", "anóm%", "r_tod", "r_limp", "disp%", "razão", "rmse%"))
for r in resultados:
    f = lambda k, fmt: (fmt % r[k]) if r.get(k) is not None else "—"
    print("%-11s %8.0f %6.1f %5s %6s %6s %6s %7s %7s %7s" % (
        r["nome"], r["envolvente_ha"], r["margem_km_por_km2"], r.get("datas_utilizaveis", 0),
        f("pct_anomalas", "%.0f"), f("r_todas", "%.2f"), f("r_sem_anomalas", "%.2f"),
        f("dispersao_pct", "%.1f"), f("razao_area_satelite_oficial", "%.2f"),
        f("serie_rmse_pct_volume_medio", "%.1f")))

# --- teste da hipótese da margem --------------------------------------------
from math import log10

excl = [r["nome"] for r in resultados if r.get("curva_oficial_consistente") is False]
par = [(r["margem_km_por_km2"], r["razao_area_satelite_oficial"], r["envolvente_ha"], r["nome"])
       for r in resultados
       if r.get("razao_area_satelite_oficial") is not None and r.get("curva_oficial_consistente")]
teste = {}
if len(par) >= 4:
    xs, ys = [p[0] for p in par], [p[1] for p in par]
    b = polyfit(xs, ys, 1)
    teste = {"n_albufeiras": len(par), "excluidas_curva_oficial_inconsistente": excl,
             "r_margem_vs_razao": round(pearson(xs, ys), 3),
             "r_log_area_vs_razao": round(pearson([log10(p[2]) for p in par], ys), 3),
             "declive_por_km_km2": round(b[1], 4), "ordenada": round(b[0], 3),
             "largura_equivalente_m": {p[3]: round(1000 * (1 - p[1]) / p[0]) for p in par}}
    print("\nHIPÓTESE DA MARGEM (%d albufeiras; excluídas por curva oficial inconsistente: %s):"
          % (len(par), ", ".join(excl) or "nenhuma"))
    print("  r (margem por área vs razão satélite/oficial) = %.2f   (hipótese: negativo)" % teste["r_margem_vs_razao"])
    print("  r (log da área vs razão)                      = %.2f" % teste["r_log_area_vs_razao"])
    print("  largura de margem que explicaria o défice (m): " +
          ", ".join("%s %d" % kv for kv in teste["largura_equivalente_m"].items()))
    print("  (um pixel tem 20 m: larguras muito acima disso não se explicam só com pixels mistos)")

campos = ["id", "nome", "envolvente_ha", "margem_km_por_km2", "datas_total", "datas_utilizaveis",
          "datas_anomalas", "pct_anomalas", "r_todas", "r_sem_anomalas", "r2_quadratico",
          "dispersao_ha", "dispersao_pct", "razao_area_satelite_oficial", "razao_dv_satelite_oficial",
          "curva_oficial_consistente", "serie_rmse_hm3",
          "serie_rmse_pct_volume_medio", "serie_rmse_pct_variacao", "nota"]
with open(os.path.join("dados", "resumo_fase3.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=campos, extrasaction="ignore")
    w.writeheader()
    w.writerows(resultados)
with open(os.path.join("dados", "resumo_fase3.json"), "w", encoding="utf-8") as fh:
    json.dump({"albufeiras": resultados, "hipotese_margem": teste,
               "fontes": "Área: Sentinel-2 (análise própria). Cota e volume: SNIRH (APA, I.P.). "
                         "Envolventes: JRC Global Surface Water."}, fh, ensure_ascii=False, indent=2)
print("\nGuardado: dados/resumo_fase3.csv, dados/resumo_fase3.json e, por albufeira, area_cota_ajuste.csv, volume_comparacao.csv, resumo.json")
