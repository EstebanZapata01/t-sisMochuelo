#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tabla de A_90(Ar, N_e en 1-5, 62 kg dia) frente al NIVEL de fondo declarado
por la ref.[46] (B_k de bkg_bins_Ar.dat, flat por bin) y a su normalizacion:
misma formula de Asimov perfilada que fig_fondo_perfilado.py (nuisance
gaussiano lineal beta, B_k(beta)=(1+beta) B_k^0),
    S2 = sum_k S_k^2/(S_k+B_k),  P = sum_k S_k B_k/(S_k+B_k),
    Q = sum_k B_k^2/(S_k+B_k),   S2_eff(sigma_B) = S2 - P^2/(Q+1/sigma_B^2),
    A_90 = 1 + sqrt(2.706/S2_eff),
evaluada en tres normalizaciones (sigma_B=0: fondo conocido; sigma_B=1%;
sigma_B->infinito: fondo libre) para cuatro niveles de fondo (0, 0.1, 1 y
10 veces el nivel declarado). Nivel x0 = sin fondo (techo intrinseco); no
depende de sigma_B, se reporta en una sola columna.

Entrada: datos/bkg_bins_Ar.dat (N_e, S_k, B_k, ev/(kg dia); chi2_bkg_nest.f90)
Salida : stdout + datos/tabla_fondo_Ar.tex
"""
import numpy as np

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
E0 = 62.0          # kg dia, exposicion de la comparacion (bkg_bins_Ar.dat)
DCHI2_90 = 2.706    # 1 g.d.l., intervalo de dos lados al 90% C.L. (el unilateral seria 1.642)

Ne, Sk, Bk0 = np.loadtxt(f"{BASE}/bkg_bins_Ar.dat", comments="#", unpack=True)
S = E0 * Sk

def A90(S2eff):
    return 1.0 + np.sqrt(DCHI2_90 / S2eff)

NIVELES = [0.0, 0.1, 1.0, 10.0]
SIGMA_B_1PCT = 0.01

filas = []
print("=" * 60)
print(f"{'fondo':>8s} {'conocida':>10s} {'+-1%':>10s} {'libre':>10s}")
print("-" * 60)
for mult in NIVELES:
    if mult == 0.0:
        S2 = np.sum(S**2 / S)
        a = A90(S2)
        filas.append((mult, a, None, None))
        print(f"{'0 (sin)':>8s} {a:10.4f} {'--':>10s} {'--':>10s}")
        continue
    B = E0 * Bk0 * mult
    sig2 = S + B
    S2 = np.sum(S**2 / sig2)
    P = np.sum(S * B / sig2)
    Q = np.sum(B**2 / sig2)
    a_conocida = A90(S2)
    a_1pct = A90(S2 - P**2 / (Q + 1.0 / SIGMA_B_1PCT**2))
    a_libre = A90(S2 - P**2 / Q)
    filas.append((mult, a_conocida, a_1pct, a_libre))
    print(f"{mult:7.1f}x {a_conocida:10.4f} {a_1pct:10.4f} {a_libre:10.4f}")
print("-" * 60)

# ------------------------------------------------------------------ .tex
with open(f"{BASE}/tabla_fondo_Ar.tex", "w", encoding="utf-8") as f:
    f.write("% generado por python/tabla_fondo_Ar.py\n")
    f.write("\\begin{tabular}{@{}lccc@{}}\n\\toprule\n")
    f.write("Nivel de fondo & conocida ($\\sigma_B{=}0$) & $\\sigma_B{=}1\\%$ & "
            "libre ($\\sigma_B\\to\\infty$) \\\\\n\\midrule\n")
    for mult, a_c, a_1, a_l in filas:
        if mult == 0.0:
            f.write(f"sin fondo & \\multicolumn{{3}}{{c}}{{${a_c:.3f}$}} \\\\\n")
        else:
            f.write(f"${mult:g}\\times$ declarado & ${a_c:.3f}$ & ${a_1:.3f}$ & ${a_l:.3f}$ \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print(f"\n  {BASE}/tabla_fondo_Ar.tex")
