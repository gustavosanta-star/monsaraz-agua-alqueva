"""
Monsaraz em Mudança — passo 10
DADOS OFICIAIS DO SNIRH (APA): cota e volume armazenado da albufeira do Alqueva.

Estes dados são medidos no terreno (estação 24M/07A), por isso são uma
referência INDEPENDENTE das imagens de satélite. São a base da validação
área–cota (passo 11) e do cálculo de volume dV = A · dh (fases seguintes).

Fonte: SNIRH — Sistema Nacional de Informação de Recursos Hídricos (APA, I.P.),
https://snirh.apambiente.pt. Ao usar estes dados, referir a origem.

Produz:
  - snirh_alqueva_bruto.csv   todos os registos (data, hora, parâmetro, valor, flag)
  - snirh_alqueva_diario.csv  cota média diária + n.º de registos; volume do dia
  - snirh_alqueva_t3.csv      resumo da cota no 3.º trimestre de cada ano
"""

import csv
import re
import urllib.request
from collections import defaultdict
from statistics import mean

BASE = "https://snirh.apambiente.pt/snirh/_dadosbase/site/paraCSV/dados_csv.php"

# Estação e parâmetros (códigos internos do SNIRH)
ESTACAO = {"site": "1627743416", "codigo": "24M/07A", "nome": "Albufeira do Alqueva (R.E.)"}
PARAMETROS = {
    "cota_m": "1629599726",      # Cota da albufeira (m)
    "volume_dam3": "1866",       # Volume armazenado (convencional), em dam³ (1000 m³)
}
TMIN, TMAX = "01/01/2017", "31/12/2025"

# linha de dados: 31/01/2017 23:00,147.57,(vco),
LINHA = re.compile(r"^(\d{2})/(\d{2})/(\d{4}) (\d{2}:\d{2}),([-\d.]+),([^,]*),")


def descarregar(par_id):
    url = (f"{BASE}?sites={ESTACAO['site']}&pars={par_id}"
           f"&tmin={TMIN}&tmax={TMAX}&formato=csv")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (monsaraz-agua-alqueva)"})
    with urllib.request.urlopen(req, timeout=180) as r:
        texto = r.read().decode("latin-1")          # o SNIRH serve em latin-1
    registos = []
    for linha in texto.splitlines():
        m = LINHA.match(linha)
        if m:
            d, mes, ano, hora, valor, flag = m.groups()
            registos.append((f"{ano}-{mes}-{d}", hora, float(valor), flag.strip().strip("()")))
    return registos


print("A descarregar do SNIRH: %s (%s), %s a %s" % (ESTACAO["nome"], ESTACAO["codigo"], TMIN, TMAX))
dados = {}
for nome, par_id in PARAMETROS.items():
    dados[nome] = descarregar(par_id)
    print("  %-12s %6d registos" % (nome, len(dados[nome])))

# --- 1. bruto ---------------------------------------------------------------
with open("snirh_alqueva_bruto.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["data", "hora", "parametro", "valor", "flag"])
    for nome, regs in dados.items():
        for data, hora, valor, flag in regs:
            w.writerow([data, hora, nome, valor, flag])

# --- 2. diário --------------------------------------------------------------
cota_dia = defaultdict(list)
for data, _, valor, _ in dados["cota_m"]:
    cota_dia[data].append(valor)
vol_dia = {}
for data, _, valor, _ in dados["volume_dam3"]:
    vol_dia[data] = valor                            # último registo do dia

datas = sorted(set(cota_dia) | set(vol_dia))
with open("snirh_alqueva_diario.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["data", "cota_m", "n_registos_cota", "volume_dam3"])
    for d in datas:
        cs = cota_dia.get(d)
        w.writerow([d,
                    round(mean(cs), 3) if cs else "",
                    len(cs) if cs else 0,
                    vol_dia.get(d, "")])

# --- 3. resumo do T3 (julho–setembro) ---------------------------------------
print("\nCota no 3.º trimestre (julho–setembro):")
print("ano | média (m) | mín (m) | máx (m) | dias c/ dados")
print("----+-----------+---------+---------+--------------")
linhas = []
for ano in range(int(TMIN[-4:]), int(TMAX[-4:]) + 1):
    dias = [d for d in cota_dia if d.startswith(str(ano)) and 7 <= int(d[5:7]) <= 9]
    if not dias:
        print("%4d| sem dados" % ano)
        continue
    medias = [mean(cota_dia[d]) for d in dias]
    linhas.append([ano, round(mean(medias), 2), round(min(medias), 2), round(max(medias), 2), len(dias)])
    aviso = "  <-- poucos dias" if len(dias) < 30 else ""
    print("%4d| %9.2f | %7.2f | %7.2f | %12d%s" % (ano, linhas[-1][1], linhas[-1][2],
                                                   linhas[-1][3], len(dias), aviso))

with open("snirh_alqueva_t3.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["ano", "cota_media_m", "cota_min_m", "cota_max_m", "dias_com_dados"])
    w.writerows(linhas)

print("\nGuardado: snirh_alqueva_bruto.csv, snirh_alqueva_diario.csv, snirh_alqueva_t3.csv")
print("Fonte: SNIRH (APA, I.P.) — https://snirh.apambiente.pt")
