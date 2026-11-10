"""
Monsaraz em Mudança — passo 12 (Fases 1 e 2)
VALIDAÇÃO ÁREA–COTA E VOLUME, a partir dos ficheiros já produzidos.
Não usa o Earth Engine: corre em segundos e só precisa de Python.

Fase 1 — curva área–cota (validação independente):
  - lê area_cota_alqueva.csv (passo 11): área do Alqueva por data + cota oficial;
  - usa as datas com o Alqueva quase todo visível (lacuna < 50 ha) e cota medida;
  - ajusta A(h) = curva de 2.º grau, com rejeição robusta de datas anómalas
    (desvio > 3,5 × dispersão robusta, MAD), repetida até estabilizar;
  - as datas rejeitadas são prováveis erros de classificação (neblina, nuvens finas,
    água turva) e ficam assinaladas, não apagadas.

Fase 2 — volume:
  - lê snirh_alqueva_diario.csv (passo 10): cota e volume oficiais;
  - ajusta V(h) oficial (3.º grau), com a mesma rejeição robusta de registos anómalos;
  - compara a área do satélite com a área implícita no volume oficial (dV/dh);
  - estima a variação de volume por satélite: ΔV = ∫ A(h) dh;
  - reconstrói a série de volume a partir de UM volume de referência
    (o satélite mede variações de volume, não o volume total: não vê o fundo).

Produz:
  - area_cota_ajuste.csv     por data: cota, área, área prevista, resíduo, anómala (sim/não)
  - volume_comparacao.csv    por registo oficial: cota, volume oficial, volume estimado
  - resumo_validacao.json    números-chave (para a app e para o README)
"""

import csv
import json
from statistics import mean, median, pstdev

LACUNA_MAX_HA = 50
K_MAD = 3.5


# --- utilitários (sem numpy) -------------------------------------------------
def polyfit(xs, ys, grau):
    """Mínimos quadrados: coeficientes [c0, c1, ...] de c0 + c1·x + c2·x² + ..."""
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
    """Integral de a a b do polinómio c."""
    return sum(ci / (i + 1) * (b ** (i + 1) - a ** (i + 1)) for i, ci in enumerate(c))


def polyder(c):
    return [i * ci for i, ci in enumerate(c)][1:]


def ajuste_robusto(xs, ys, grau, x0):
    """Ajuste polinomial em (x - x0) com rejeição iterativa de anómalos (MAD)."""
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


# --- Fase 1: curva área–cota -------------------------------------------------
with open("area_cota_alqueva.csv", encoding="utf-8") as fh:
    todas = list(csv.DictReader(fh))
uteis = [r for r in todas if r["cota_m"] and float(r["lacuna_ha"]) < LACUNA_MAX_HA]
H = [float(r["cota_m"]) for r in uteis]
A = [float(r["agua_ha"]) for r in uteis]
HA0 = round(mean(H), 2)
cA, resA, okA = ajuste_robusto(H, A, 2, HA0)

Hc = [h for h, k in zip(H, okA) if k]
Ac = [a for a, k in zip(A, okA) if k]
rA_todas = pearson(H, A)
rA_limpo = pearson(Hc, Ac)
res_limpo = [r for r, k in zip(resA, okA) if k]
ss_res = sum(r * r for r in res_limpo)
ss_tot = sum((a - mean(Ac)) ** 2 for a in Ac)
disp_ha = pstdev(res_limpo)

with open("area_cota_ajuste.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["data", "cota_m", "agua_ha", "prevista_ha", "residuo_ha", "anomala"])
    for r, h, a, res, k in zip(uteis, H, A, resA, okA):
        w.writerow([r["data"], h, a, round(a - res, 1), round(res, 1), "nao" if k else "sim"])

print("FASE 1 — curva área–cota do Alqueva")
print("  datas utilizáveis: %d (lacuna < %d ha e cota medida)" % (len(H), LACUNA_MAX_HA))
print("  na curva: %d | anómalas: %d (%.0f%%)" % (len(Hc), len(H) - len(Hc), 100 * (len(H) - len(Hc)) / len(H)))
print("  r com todas as datas: %.3f" % rA_todas)
print("  r sem as anómalas:    %.3f | R² (2.º grau) = %.3f" % (rA_limpo, 1 - ss_res / ss_tot))
print("  dispersão: ±%.0f ha (%.1f%% da área média)" % (disp_ha, 100 * disp_ha / mean(Ac)))

# --- Fase 2: volume ----------------------------------------------------------
with open("snirh_alqueva_diario.csv", encoding="utf-8") as fh:
    of = [r for r in csv.DictReader(fh) if r["cota_m"] and r["volume_dam3"]]
Ho = [float(r["cota_m"]) for r in of]
Vo = [float(r["volume_dam3"]) / 1000 for r in of]           # dam³ → hm³
HV0 = 147.5
cV, resV, okV = ajuste_robusto(Ho, Vo, 3, HV0)
dispV = pstdev([r for r, k in zip(resV, okV) if k])
anom_of = [(r["data"], h, v) for r, h, v, k in zip(of, Ho, Vo, okV) if not k]

print("\nFASE 2 — volume")
print("  registos oficiais com cota e volume: %d | anómalos: %d" % (len(Ho), len(anom_of)))
for d, h, v in anom_of:
    print("    %s: cota %.2f m, volume %.0f hm³ (assinalado)" % (d, h, v))
print("  dispersão do volume oficial em torno de V(h): %.0f hm³" % dispV)

# área do satélite vs área implícita no volume oficial (dV/dh, em ha: hm³/m × 100)
dV = polyder(cV)
comparacao = []
print("\n  cota  | área satélite | área implícita (dV/dh oficial) | satélite/oficial")
for h in [144, 145, 146, 147, 148, 149, 150, 151]:
    a_sat = polyval(cA, h - HA0)
    a_of = polyval(dV, h - HV0) * 100
    comparacao.append({"cota_m": h, "area_satelite_ha": round(a_sat), "area_oficial_ha": round(a_of),
                       "razao": round(a_sat / a_of, 3)})
    print("  %5.1f | %9.0f ha | %12.0f ha | %.2f" % (h, a_sat, a_of, a_sat / a_of))
razao_media = mean(c["razao"] for c in comparacao)

# variação de volume entre o mínimo e o máximo do período
Hok = [h for h, k in zip(Ho, okV) if k]
h_min, h_max = min(Hok), max(Hok)
dv_sat = polyint(cA, h_min - HA0, h_max - HA0) / 100     # ha·m → hm³
dv_of = polyval(cV, h_max - HV0) - polyval(cV, h_min - HV0)
print("\n  de %.2f m a %.2f m: ΔV satélite = %.0f hm³ | ΔV oficial = %.0f hm³ (razão %.2f)"
      % (h_min, h_max, dv_sat, dv_of, dv_sat / dv_of))

# série de volume reconstruída a partir de um volume de referência (registo mais perto de 147,5 m)
i_ref = min((i for i in range(len(Ho)) if okV[i]), key=lambda i: abs(Ho[i] - 147.5))
v_sat = [Vo[i_ref] + polyint(cA, Ho[i_ref] - HA0, h - HA0) / 100 for h in Ho]
erros = [vs - vo for vs, vo, k in zip(v_sat, Vo, okV) if k]
rmse = (sum(e * e for e in erros) / len(erros)) ** 0.5
v_medio = mean(v for v, k in zip(Vo, okV) if k)
print("  série (%d registos, referência %s): erro médio %+.0f hm³ | erro quadrático médio %.0f hm³ (%.1f%%)"
      % (len(erros), of[i_ref]["data"], mean(erros), rmse, 100 * rmse / v_medio))

with open("volume_comparacao.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["data", "cota_m", "volume_oficial_hm3", "volume_satelite_hm3", "anomalo_oficial"])
    for r, h, vo, vs, k in zip(of, Ho, Vo, v_sat, okV):
        w.writerow([r["data"], h, round(vo, 1), round(vs, 1), "nao" if k else "sim"])

# --- resumo para a app e o README --------------------------------------------
resumo = {
    "fase1_area_cota": {
        "datas_utilizaveis": len(H), "datas_na_curva": len(Hc), "datas_anomalas": len(H) - len(Hc),
        "r_todas": round(rA_todas, 3), "r_sem_anomalas": round(rA_limpo, 3),
        "r2_quadratico": round(1 - ss_res / ss_tot, 3),
        "dispersao_ha": round(disp_ha), "dispersao_pct": round(100 * disp_ha / mean(Ac), 1),
        "cota_min_m": min(H), "cota_max_m": max(H),
        "curva": {"forma": "A = c0 + c1*(h-h0) + c2*(h-h0)^2", "h0": HA0, "c": [round(x, 3) for x in cA]},
    },
    "fase2_volume": {
        "registos_oficiais": len(Ho), "registos_anomalos": len(anom_of),
        "razao_area_satelite_oficial_media": round(razao_media, 3), "comparacao_por_cota": comparacao,
        "dv_satelite_hm3": round(dv_sat), "dv_oficial_hm3": round(dv_of),
        "cota_min_m": h_min, "cota_max_m": h_max,
        "serie_rmse_hm3": round(rmse), "serie_rmse_pct": round(100 * rmse / v_medio, 1),
        "serie_erro_medio_hm3": round(mean(erros)), "referencia": of[i_ref]["data"],
    },
    "fontes": "Área: Sentinel-2 (análise própria). Cota e volume: SNIRH (APA, I.P.).",
}
with open("resumo_validacao.json", "w", encoding="utf-8") as fh:
    json.dump(resumo, fh, ensure_ascii=False, indent=2)
print("\nGuardado: area_cota_ajuste.csv, volume_comparacao.csv, resumo_validacao.json")
