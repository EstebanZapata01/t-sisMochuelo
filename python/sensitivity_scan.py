#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analisis de sensibilidad: cuanto se mueve el resultado ideal (A_90 y la
ventaja Ar/Xe) al variar UN parametro a la vez dentro de su incertidumbre
publicada / justificable, y si la conclusion central ("Ar gana con umbral
alto") sobrevive a todas las variaciones.

Corre chi2_ideal (Xe y Ar) bajo cada variacion via variables de entorno
(IDEAL_QY_SCALE, IDEAL_TMAX_KEV, USE_HELM); la fluctuacion de ionizacion
F (Fano vs Poisson) ya sale de las dos filas F_mode de sensib_ideal_*.
Restaura los ficheros baseline al final.

Parametros y su fuente:
  Q_y (yield)   Xe x{0.85,1.15}  (banda NR de NEST a ~keV)
                Ar x{0.92,1.08}  (razon ReD/LArNEST, arXiv:2510.16404)
  F ionizacion  Fano (NEST)  <->  Poisson (F=1)
  T_max (Ar)    3.44 keV (ref.[46], cinematica E_nu=8) <-> 5.4 keV (E_nu=10)
  Helm F^2      F^2=1 (base)  <->  Helm activado (USE_HELM=1)

Salidas: datos/sensib_tornado.dat, datos/fig_sensib_tornado.png,
         datos/tabla_sensib.tex
"""
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DIRX = f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe"
DIRA = f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"
NE_FOCUS = 4          # umbral en el que se hace el tornado
C_XE, C_AR = "#33546e", "#a86a43"


def run(target_dir, env_extra):
    e = dict(os.environ)
    e.update(env_extra)
    subprocess.run(["./chi2_ideal"], cwd=target_dir, env=e,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


def read_sensib(tag):
    """-> {(NE_LO, F_mode): (Ntot@192, A_90_x1)}"""
    out = {}
    with open(f"{BASE}/sensib_ideal_{tag}.dat", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            out[(int(p[0]), p[1])] = (float(p[4]), float(p[5]))
    return out


def a90(tag, ne, fmode="Fnest"):
    return read_sensib(tag)[(ne, fmode)][1]


def scenario(name, envx, enva, fmode_xe="Fnest", fmode_ar="Fnest"):
    """Corre Xe y Ar con el entorno dado y devuelve (A90_Xe, A90_Ar, ratio_Rtot)."""
    run(DIRX, envx)
    run(DIRA, enva)
    sx, sa = read_sensib("Xe"), read_sensib("Ar")
    ax = sx[(NE_FOCUS, fmode_xe)][1]
    aa = sa[(NE_FOCUS, fmode_ar)][1]
    # ventaja = cociente de N_tot (proporcional a R_tot) en el umbral foco
    ratio = sa[(NE_FOCUS, fmode_ar)][0] / sx[(NE_FOCUS, fmode_xe)][0]
    return ax, aa, ratio


# ----------------------------------------------------------------- baseline
run(DIRX, {}); run(DIRA, {})
base_ax, base_aa, base_ratio = scenario("baseline", {}, {})
print("=" * 82)
print(f" ANALISIS DE SENSIBILIDAD  (umbral foco: N_e >= {NE_FOCUS}, comparacion ideal)")
print("=" * 82)
print(f" baseline: A_90 Xe = {base_ax:.4f}   A_90 Ar = {base_aa:.4f}   "
      f"ventaja Ar/Xe = x{base_ratio:.1f}")
print("-" * 82)

VARIATIONS = [
    ("Q_y yield (Xe ±15%, Ar ±8%)",
     dict(lo=(dict(IDEAL_QY_SCALE="0.85"), dict(IDEAL_QY_SCALE="0.92")),
          hi=(dict(IDEAL_QY_SCALE="1.15"), dict(IDEAL_QY_SCALE="1.08")))),
    ("T_max Ar (3.44 ↔ 5.4 keV)",
     dict(lo=({}, dict(IDEAL_TMAX_KEV="3.44")),
          hi=({}, dict(IDEAL_TMAX_KEV="5.4")))),
    ("Helm F² (1 ↔ Helm)",
     dict(lo=({}, {}),
          hi=(dict(USE_HELM="1"), dict(USE_HELM="1")))),
]

rows = []   # (label, A90_Xe_lo, A90_Xe_hi, A90_Ar_lo, A90_Ar_hi, ratio_lo, ratio_hi)
for label, vv in VARIATIONS:
    xl, al, rl = scenario(label + " lo", *vv["lo"])
    xh, ah, rh = scenario(label + " hi", *vv["hi"])
    rows.append((label, xl, xh, al, ah, rl, rh))
    print(f" {label:32s}  A90_Xe [{xl:.3f}, {xh:.3f}]  A90_Ar [{al:.4f}, {ah:.4f}]  "
          f"Ar/Xe [{rl:.0f}, {rh:.0f}]")

# F Fano <-> Poisson desde el baseline (ambas filas ya calculadas)
run(DIRX, {}); run(DIRA, {})
sx, sa = read_sensib("Xe"), read_sensib("Ar")
f_ax = (sx[(NE_FOCUS, "Fnest")][1], sx[(NE_FOCUS, "F1")][1])
f_aa = (sa[(NE_FOCUS, "Fnest")][1], sa[(NE_FOCUS, "F1")][1])
f_r  = (sa[(NE_FOCUS, "Fnest")][0] / sx[(NE_FOCUS, "Fnest")][0],
        sa[(NE_FOCUS, "F1")][0]    / sx[(NE_FOCUS, "F1")][0])
rows.append(("F ionización (Fano ↔ Poisson)",
             min(f_ax), max(f_ax), min(f_aa), max(f_aa), min(f_r), max(f_r)))
print(f" {'F ionización (Fano ↔ Poisson)':32s}  A90_Xe [{min(f_ax):.3f}, {max(f_ax):.3f}]  "
      f"A90_Ar [{min(f_aa):.4f}, {max(f_aa):.4f}]  Ar/Xe [{min(f_r):.0f}, {max(f_r):.0f}]")

# ----------------------------------------------------------------- restaurar
run(DIRX, {}); run(DIRA, {})
print("-" * 82)
print(" baseline restaurado.")

# ----------------------------------------------------------------- tornado fig
labels = [r[0] for r in rows]
y = np.arange(len(rows))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 0.75 * len(rows) + 2.6))

for ax, idx_lo, idx_hi, cbase, title in (
        (a1, 1, 2, C_XE, r"$A_{90}$ de Xe ($N_e\geq4$, ideal)"),
        (a2, 5, 6, "0.6", r"ventaja Ar/Xe ($N_e\geq4$)")):
    base = base_ax if ax is a1 else base_ratio
    for i, r in enumerate(rows):
        lo, hi = sorted((r[idx_lo], r[idx_hi]))
        ax.barh(i, hi - lo, left=lo, height=0.55,
                color=cbase, alpha=0.55, edgecolor="black", lw=0.7)
        ax.plot([lo, hi], [i, i], "|", color="black", ms=10, mew=1.3)
    ax.axvline(base, color="black", lw=1.3, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels if ax is a1 else [])
    ax.set_title(title, pad=16)
    ax.annotate(f"base = {base:.3g}", xy=(base, -0.7), ha="center", va="bottom",
                fontsize=8.5, annotation_clip=False)
    ax.grid(axis="x", ls=":", lw=0.4, color="0.8")
    ax.set_ylim(len(rows) - 0.4, -1.1)

fig.suptitle("Análisis de sensibilidad: rango de $A_{90}$ y de la ventaja Ar/Xe "
             "al variar cada parámetro (uno a la vez)", y=1.0, fontsize=12)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_sensib_tornado.png")

with open(f"{BASE}/sensib_tornado.dat", "w", encoding="utf-8") as f:
    f.write(f"# umbral N_e>={NE_FOCUS} | baseline A90_Xe={base_ax:.5f} "
            f"A90_Ar={base_aa:.5f} ratio={base_ratio:.4f}\n")
    f.write("# parametro | A90_Xe_lo A90_Xe_hi | A90_Ar_lo A90_Ar_hi | ratio_lo ratio_hi\n")
    for r in rows:
        f.write(f"{r[0]:34s} {r[1]:9.5f} {r[2]:9.5f}  {r[3]:9.5f} {r[4]:9.5f}  "
                f"{r[5]:9.3f} {r[6]:9.3f}\n")

# etiquetas seguras para LaTeX (las de arriba llevan _, %, unicode)
TEX_LABEL = {
    "Q_y yield (Xe ±15%, Ar ±8%)":
        r"$Q_y$ (Xe $\pm15\%$, Ar $\pm8\%$)",
    "T_max Ar (3.44 ↔ 5.4 keV)":
        r"$T_{\max}$ Ar ($3{,}44\leftrightarrow5{,}4$ keV)",
    "Helm F² (1 ↔ Helm)":
        r"Factor de forma ($F^2{=}1\leftrightarrow$ Helm)",
    "F ionización (Fano ↔ Poisson)":
        r"Fluctuaci\'on $F$ (Fano $\leftrightarrow$ Poisson)",
}
with open(f"{BASE}/tabla_sensib.tex", "w", encoding="utf-8") as f:
    f.write(f"% python/sensitivity_scan.py  (umbral N_e = {NE_FOCUS})\n")
    f.write("\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{5.4cm}ccc@{}}"
            "\n\\toprule\n")
    f.write("Par\\'ametro variado & rango $A_{90}$ Xe & rango $A_{90}$ Ar & "
            "rango Ar/Xe \\\\\n\\midrule\n")
    f.write(f"\\emph{{baseline}} & {base_ax:.3f} & {base_aa:.4f} & "
            f"${base_ratio:.0f}\\times$ \\\\\n")
    for r in rows:
        lab = TEX_LABEL.get(r[0], r[0])
        f.write(f"{lab} & $[{min(r[1],r[2]):.3f},{max(r[1],r[2]):.3f}]$ & "
                f"$[{min(r[3],r[4]):.4f},{max(r[3],r[4]):.4f}]$ & "
                f"$[{min(r[5],r[6]):.0f},{max(r[5],r[6]):.0f}]\\times$ \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# conclusion robustez
all_ratio_lo = min(r[5] for r in rows)
print(f"\n  Peor caso de la ventaja Ar/Xe (N_e>={NE_FOCUS}) sobre todas las "
      f"variaciones: x{all_ratio_lo:.0f}  (baseline x{base_ratio:.0f})")
print(f"  -> la conclusión 'Ar gana con umbral alto' "
      f"{'SOBREVIVE' if all_ratio_lo > 5 else 'NO sobrevive'} a todas las variaciones.")
print(f"\n  {BASE}/fig_sensib_tornado.png\n  {BASE}/sensib_tornado.dat\n  {BASE}/tabla_sensib.tex")
