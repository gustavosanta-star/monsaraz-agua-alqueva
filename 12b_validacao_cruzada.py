"""
Monsaraz em Mudança — passo 12b
TESTE MAIS EXIGENTE DO VOLUME: satélite sozinho e fora da amostra.

O passo 12 reconstrói o volume com a curva A(h) do satélite mas usando a COTA OFICIAL
de cada dia, e avalia a curva nas mesmas datas em que foi ajustada. Aqui testa-se o uso
real: depois de calibrada, UMA imagem dá o volume sem qualquer medição no terreno.

  área da imagem  ->  cota estimada (curva A(h) invertida)  ->  volume (V(h) oficial)

e compara-se com o volume oficial desse dia. Dois testes:
  1. na amostra (curva ajustada com todas as datas);
  2. deixando um ano de fora: a curva é ajustada sem esse ano e testada nele
     (as datas anómalas do ano de teste são definidas pela curva dos outros anos).

Não usa o Earth Engine. Requer area_cota_alqueva.csv (passo 11) e snirh_alqueva_diario.csv (passo 10).
Produz: validacao_cruzada.json
"""

import csv
import json
import random
from statistics import mean, median

LACUNA_MAX_HA = 50
K_MAD = 3.5


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
    return c, med, mad, ok


def pearson(xs, ys):
    mx, my = mean(xs), mean(ys)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    return sxy / (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5


def rmse(e):
    return (sum(x * x for x in e) / len(e)) ** 0.5


with open("area_cota_alqueva.csv", encoding="utf-8") as fh:
    rows = [r for r in csv.DictReader(fh) if r["cota_m"] and float(r["lacuna_ha"]) < LACUNA_MAX_HA]
D = [r["data"] for r in rows]
H = [float(r["cota_m"]) for r in rows]
A = [float(r["agua_ha"]) for r in rows]
H0 = round(mean(H), 2)

with open("snirh_alqueva_diario.csv", encoding="utf-8") as fh:
    of = [r for r in csv.DictReader(fh) if r["cota_m"] and r["volume_dam3"]]
HV0 = 147.5
cV = ajuste_robusto([float(r["cota_m"]) for r in of], [float(r["volume_dam3"]) / 1000 for r in of], 3, HV0)[0]
V = lambda h: polyval(cV, h - HV0)                      # volume oficial (hm³) à cota h


def cota_da_area(c, a, lo=138.0, hi=156.0):
    """Inverte a curva A(h) (crescente) por bissecção."""
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if polyval(c, m - H0) < a else (lo, m)
    return (lo + hi) / 2


def erros(c, idx):
    eh = [cota_da_area(c, A[i]) - H[i] for i in idx]
    ev = [V(cota_da_area(c, A[i])) - V(H[i]) for i in idx]
    return eh, ev


def resumo(eh, ev, idx):
    vm = mean(V(H[i]) for i in idx)
    return {"n": len(idx), "cota_vies_m": round(mean(eh), 2), "cota_rmse_m": round(rmse(eh), 2),
            "volume_vies_hm3": round(mean(ev)), "volume_rmse_hm3": round(rmse(ev)),
            "volume_rmse_pct": round(100 * rmse(ev) / vm, 1),
            "volume_erro_mediano_hm3": round(median(abs(x) for x in ev))}


todos = list(range(len(H)))
cA, med, mad, ok = ajuste_robusto(H, A, 2, H0)
bons = [i for i in todos if ok[i]]
na_amostra = resumo(*erros(cA, bons), bons)
com_anomalas = resumo(*erros(cA, todos), todos)

por_ano, eh_f, ev_f, idx_f = [], [], [], []
for ano in sorted(set(d[:4] for d in D)):
    tr = [i for i in todos if D[i][:4] != ano]
    te = [i for i in todos if D[i][:4] == ano]
    c, m, s, _ = ajuste_robusto([H[i] for i in tr], [A[i] for i in tr], 2, H0)
    te_ok = [i for i in te if abs(A[i] - polyval(c, H[i] - H0) - m) <= K_MAD * s]
    eh, ev = erros(c, te_ok)
    por_ano.append({"ano": int(ano), "datas": len(te), "na_curva": len(te_ok),
                    "volume_vies_hm3": round(mean(ev)), "volume_rmse_hm3": round(rmse(ev))})
    eh_f += eh; ev_f += ev; idx_f += te_ok
fora = resumo(eh_f, ev_f, idx_f)

random.seed(1)
rs = sorted(pearson([H[i] for i in s], [A[i] for i in s])
            for s in ([random.choice(bons) for _ in bons] for _ in range(2000)))
ic = [round(rs[50], 3), round(rs[1949], 3)]

print("VOLUME A PARTIR DE UMA SÓ IMAGEM (área -> cota -> volume), Alqueva")
print("  na amostra, sem anómalas : n=%(n)d | cota ±%(cota_rmse_m).2f m | volume ±%(volume_rmse_hm3)d hm³ (%(volume_rmse_pct).1f%%)" % na_amostra)
print("  FORA da amostra (1 ano de fora): n=%(n)d | cota ±%(cota_rmse_m).2f m | volume ±%(volume_rmse_hm3)d hm³ (%(volume_rmse_pct).1f%%), viés %(volume_vies_hm3)+d hm³" % fora)
print("  na amostra, COM anómalas : n=%(n)d | volume ±%(volume_rmse_hm3)d hm³ (%(volume_rmse_pct).1f%%)" % com_anomalas)
print("\n  ano  | datas | na curva | viés (hm³) | erro (hm³)")
for p in por_ano:
    print("  %(ano)d | %(datas)5d | %(na_curva)8d | %(volume_vies_hm3)+10d | %(volume_rmse_hm3)10d" % p)
print("\n  r (área vs cota, sem anómalas): intervalo de 95%% por bootstrap = %.3f–%.3f" % tuple(ic))

with open("validacao_cruzada.json", "w", encoding="utf-8") as fh:
    json.dump({"descricao": "Volume estimado por uma só imagem (área -> cota pela curva A(h) invertida -> volume pela V(h) oficial), comparado com o volume oficial do dia.",
               "na_amostra_sem_anomalas": na_amostra, "fora_da_amostra_um_ano_de_fora": fora,
               "na_amostra_com_anomalas": com_anomalas, "por_ano": por_ano,
               "r_area_cota_ic95_bootstrap": ic,
               "fontes": "Área: Sentinel-2 (análise própria). Cota e volume: SNIRH (APA, I.P.)."},
              fh, ensure_ascii=False, indent=2)
print("\nGuardado: validacao_cruzada.json")
