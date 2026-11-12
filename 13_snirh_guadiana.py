"""
Monsaraz em Mudança — passo 13 (Fase 3)
DADOS OFICIAIS DO SNIRH PARA TODAS AS ALBUFEIRAS DO GUADIANA (albufeiras.json).

É o passo 10 generalizado: a lista de albufeiras vem do ficheiro de configuração,
não do código. Acrescentar uma albufeira = acrescentar uma linha ao albufeiras.json.

Fonte: SNIRH — Sistema Nacional de Informação de Recursos Hídricos (APA, I.P.),
https://snirh.apambiente.pt. Ao usar estes dados, referir a origem.

Uso:
  python 13_snirh_guadiana.py            todas as albufeiras
  python 13_snirh_guadiana.py caia vigia só estas

Produz, para cada albufeira, em dados/<id>/:
  - snirh_diario.csv   data, cota_m (média diária), n_registos_cota, volume_dam3
E um resumo: dados/resumo_snirh.csv
"""

import csv
import json
import os
import re
import sys
import time
import urllib.request
from collections import defaultdict
from statistics import mean

BASE = "https://snirh.apambiente.pt/snirh/_dadosbase/site/paraCSV/dados_csv.php"
TMIN, TMAX = "01/01/2017", "31/12/2025"
LINHA = re.compile(r"^(\d{2})/(\d{2})/(\d{4}) (\d{2}:\d{2}),([-\d.]+),([^,]*),")

with open("albufeiras.json", encoding="utf-8") as fh:
    CFG = json.load(fh)
PARS = CFG["parametros_snirh"]


def descarregar(site, par_id, tentativas=3):
    url = f"{BASE}?sites={site}&pars={par_id}&tmin={TMIN}&tmax={TMAX}&formato=csv"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (monsaraz-agua-alqueva)"})
    for k in range(tentativas):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                texto = r.read().decode("latin-1")
            break
        except Exception as e:
            if k == tentativas - 1:
                raise
            print("    (nova tentativa: %s)" % str(e)[:60])
            time.sleep(5)
    regs = []
    for linha in texto.splitlines():
        m = LINHA.match(linha)
        if m:
            d, mes, ano, hora, valor, flag = m.groups()
            regs.append((f"{ano}-{mes}-{d}", float(valor)))
    return regs


pedidas = sys.argv[1:]
resumo = []
for alb in CFG["albufeiras"]:
    if pedidas and alb["id"] not in pedidas:
        continue
    pasta = os.path.join("dados", alb["id"])
    os.makedirs(pasta, exist_ok=True)
    print("%s (%s)" % (alb["nome"], alb["codigo"]))

    cota = defaultdict(list)
    for d, v in descarregar(alb["site"], PARS["cota_m"]):
        if v > 0:                                   # há registos a 0 que são falhas, não cotas
            cota[d].append(v)
    vol = {}
    if alb["volume"]:
        for d, v in descarregar(alb["site"], PARS["volume_dam3"]):
            vol[d] = v                              # último registo do dia

    datas = sorted(set(cota) | set(vol))
    with open(os.path.join(pasta, "snirh_diario.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["data", "cota_m", "n_registos_cota", "volume_dam3"])
        for d in datas:
            cs = cota.get(d)
            w.writerow([d, round(mean(cs), 3) if cs else "", len(cs) if cs else 0, vol.get(d, "")])

    medias = sorted(mean(v) for v in cota.values())
    q = lambda p: medias[int(p * (len(medias) - 1))] if medias else float("nan")
    resumo.append([alb["id"], alb["nome"], len(cota), len(vol), round(q(0.02), 2), round(q(0.98), 2),
                   round(q(0.98) - q(0.02), 2)])
    print("  dias com cota: %4d | dias com volume: %4d | cota %.2f–%.2f m (amplitude %.1f m)"
          % (len(cota), len(vol), q(0.02), q(0.98), q(0.98) - q(0.02)))

with open(os.path.join("dados", "resumo_snirh.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "nome", "dias_cota", "dias_volume", "cota_p2_m", "cota_p98_m", "amplitude_m"])
    w.writerows(resumo)
print("\nGuardado em dados/<albufeira>/snirh_diario.csv e dados/resumo_snirh.csv")
print("Fonte: SNIRH (APA, I.P.) — https://snirh.apambiente.pt")
