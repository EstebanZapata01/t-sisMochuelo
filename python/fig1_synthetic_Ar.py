#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"Fig. 1" sintetica para argon: el analogo del residuo ON-OFF real de Xe
(figs. 8/9 del documento), pero para Ar es un ASIMOV PURO -- no una
medicion. Sirve para poner Xe y Ar en el mismo formato visual:
"esto es lo que se midio en Xe; esto es lo que un run equivalente de Ar
esperaria medir".

Construccion (Asimov, sin ruido, sin semilla):
  dN_k   = S_k                      (el residuo esperado ES la senal SM)
  sigma_k = sqrt( (S_k + B_k) / T ) (estadistica de conteo a exposicion T)
con S_k (senal CEvNS SM por bin) y B_k (fondo de 39Ar depletado, UAr) que
calcula chi2_bkg_nest.f90; T = 62 kg*dia (ref.[46]).

Ajuste identico al del Xe real (forma cerrada):
  chi2(A) = S3 - 2 A S1 + A^2 S2 ,  con dN_k = S_k  =>  S1 = S2 = S3
  A_best = 1 ,  chi2_min = 0 ,  A_90 = 1 + sqrt(2.706 / S2)

Salida: datos/fig1_synthetic_Ar.png
"""
import subprocess
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()
import matplotlib.ticker as ticker

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DIRA = f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"
T_EXPO = 62.0           # kg*dia, ref.[46]
NARANJA = "#a86a43"

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.labelsize": 11,
    "axes.titlesize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 8.6, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 250, "savefig.bbox": "tight",
})

# ---------------------------------------------------------------- S_k, B_k
out = subprocess.run([f"{DIRA}/chi2_bkg"], cwd=DIRA, capture_output=True,
                     text=True, check=True).stdout
S, Buar = {}, {}
sec = None
for ln in out.splitlines():
    if "senal CEvNS SM por bin" in ln:
        sec = "S"; continue
    if "fondo 39Ar por bin" in ln:
        sec = "B"; continue
    m = re.match(r"\s*N_e\s*=\s*(\d+)\s*:\s*([-\d.E+]+)(?:\s+([-\d.E+]+))?", ln)
    if not m:
        continue
    k = int(m.group(1))
    if sec == "S":
        S[k] = float(m.group(2))
    elif sec == "B":
        Buar[k] = float(m.group(2))          # 1a col = UAr depletado
        if k == 5:
            sec = None

ks = sorted(S)
Sk = np.array([S[k] for k in ks])            # ev/(kg dia)
Bk = np.array([Buar[k] for k in ks])         # decaim/(kg dia)
sig = np.sqrt((Sk + Bk) / T_EXPO)            # error estadistico del residuo [ev/(kg dia)]
dN = Sk.copy()                               # Asimov: el residuo esperado = senal SM

S1 = np.sum(dN * Sk / sig**2)
S2 = np.sum(Sk**2 / sig**2)
S3 = np.sum(dN**2 / sig**2)
A_best = S1 / S2
chi2_min = S3 - S1**2 / S2
A_90 = A_best + np.sqrt(2.706 / S2)

print("=" * 66)
print(" FIG. 1 SINTETICA DE Ar (Asimov, UAr depletado, 62 kg*dia)")
print("=" * 66)
for k, s, b, e in zip(ks, Sk, Bk, sig):
    print(f"  N_e={k}:  S={s:8.4f}  B(UAr)={b:.3e}  sigma={e:8.4f}  ev/(kg dia)")
print(f"  S1=S2=S3 = {S2:.3f}   A_best={A_best:.4f}   chi2_min={chi2_min:.2e}")
print(f"  A_90 (Asimov) = {A_90:.5f}  x SM   (cf. sensib_bkg_Ar.dat: 1.03869)")

# ---------------------------------------------------------------- figura
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 3.9))

# --- panel A: residuo Asimov + banda de exclusion 90% C.L.
a1.axhline(0.0, color="0.55", lw=0.6, ls="--")
a1.step(ks, Sk, where="mid", color="#33546e", lw=1.5, label="Predicción CEνNS (SM)")
a1.fill_between(ks, 0, A_90 * Sk, step="mid", color=NARANJA, alpha=0.28,
               label=f"Límite 90% C.L. ($A_{{90}}={A_90:.3f}$)")
a1.step(ks, A_90 * Sk, where="mid", color=NARANJA, lw=1.4, alpha=0.8)
a1.errorbar(ks, dN, yerr=sig, fmt="ko", ms=4.5, lw=0.9, capsize=2.5,
            label="Residuo esperado (Asimov)")
a1.set_xlabel(r"$N_e$ (electrones de ionización)")
a1.set_ylabel(r"Tasa $\cdot$ kg$^{-1}\cdot$día$^{-1}$")
a1.set_title("Ar — residuo ON$-$OFF esperado (Asimov)")
a1.xaxis.set_major_locator(ticker.MultipleLocator(1))
a1.legend(loc="upper right")
a1.set_ylim(-2, max(A_90 * Sk) * 1.25)

# --- panel B: perfil Delta chi2(A)
A = np.linspace(0.90, 1.12, 600)
dchi2 = (1.0 - A) ** 2 * S2
a2.plot(A, dchi2, color="k", lw=1.8, label=r"$\Delta\chi^2(A)$")
a2.axhline(2.706, color=NARANJA, lw=1.3, ls="--", label=r"90% C.L. ($\Delta\chi^2=2{,}706$)")
a2.axvline(1.0, color="0.55", lw=0.7, ls=":")
a2.axvline(A_90, color=NARANJA, lw=0.7, ls=":")
a2.annotate(rf"$A_{{90}}={A_90:.3f}$", xy=(A_90, 2.706),
            xytext=(A_90 + 0.005, 5.0), fontsize=8.5, color=NARANJA,
            arrowprops=dict(arrowstyle="->", color=NARANJA, lw=0.7))
a2.text(1.0, 7.6, "SM\n$(A=1)$", ha="center", va="top", fontsize=8, color="0.5")
a2.set_xlabel(r"Amplitud de la señal CE$\nu$NS  $A$")
a2.set_ylabel(r"$\Delta\chi^2 = \chi^2(A) - \chi^2_{\min}$")
a2.set_title(r"Ar — perfil $\Delta\chi^2$ esperado (Asimov)")
a2.set_ylim(-0.3, 9)
a2.set_xlim(0.90, 1.12)
a2.legend(loc="upper left")

fig.suptitle("Asimov (esperado), NO medición  —  formato análogo a las figs. "
             "del Xe real, con Ar depletado (UAr) a 62 kg·día", y=1.03, fontsize=10)
fig.tight_layout()
fig.savefig(f"{BASE}/fig1_synthetic_Ar.png")
print(f"\n  {BASE}/fig1_synthetic_Ar.png")
