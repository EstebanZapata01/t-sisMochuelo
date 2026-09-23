#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONUS+ (Ge): (1) los 19 bins de exceso con la prediccion del ajuste de amplitud
y (2) el perfil Delta chi^2 de eps_ee^dV (eps_emu^dV = 0) con la solucion de
cancelacion marcada. Bin j: E = 0.16 + 0.01*(j-1/2) keV_ee (constants.f90).

Entradas: datos/generic_input_conus.dat, generic_conus_resumen.txt, chi2_nsi_2Dconus.dat
Salidas : datos/fig_conus_espectro.png, datos/fig_conus_perfil_eps.png
"""
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_GE, NEGRO, GRIS, FIG15, nota
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
b, R, sg, Rp = np.loadtxt(f"{BASE}/generic_input_conus.dat", comments="#", skiprows=2, unpack=True)
A = float(re.search(r"A_best\s*=\s*([0-9.]+)", open(f"{BASE}/generic_conus_resumen.txt").read()).group(1))
E = 160.0 + 10.0 * (b - 0.5)                       # eV_ee

fig, ax = plt.subplots(figsize=FIG15)
ax.errorbar(E, R, yerr=sg, fmt="o", color=NEGRO, ms=6, capsize=3, lw=1.2, label="exceso medido")
ax.step(np.append(E - 5, E[-1] + 5), np.append(A * Rp, A * Rp[-1]), where="post", color=C_GE,
        label="predicción del ajuste")
ax.axhline(0, color=GRIS, lw=0.7)
ax.set_xlabel(r"$E_{\rm ee}$  [eV]")
ax.set_ylabel("eventos por bin de 10 eV")
ax.legend(loc="upper right")
fig.tight_layout(); fig.savefig(f"{BASE}/fig_conus_espectro.png"); plt.close(fig)

g = np.loadtxt(f"{BASE}/chi2_nsi_2Dconus.dat", comments="#")
x, y, c = g[:, 0], g[:, 1], g[:, 2]
ux, uy = np.unique(x), np.unique(y)
fila = c.reshape(len(uy), len(ux))[np.argmin(np.abs(uy))]          # eps_emu = 0
d2 = fila - fila.min()
i2 = np.argmin(np.where(ux > 0.2, fila, np.inf))                     # segunda solucion
fig, ax = plt.subplots(figsize=FIG15)
ax.plot(ux, d2, color=C_GE)
for niv in (2.706,):
    ax.axhline(niv, color=GRIS, ls="--", lw=0.9)
ax.text(-0.98, 2.706 + 0.35, r"$90\%$ C.L.", color=GRIS, fontsize=11)
ax.plot(ux[i2], d2[i2], "o", color=NEGRO, ms=7)
nota(ax, (ux[i2], d2[i2]), "solución de\ncancelación", (ux[i2] + 0.12, d2[i2] + 4.0))
ax.set_xlim(ux[0], ux[-1]); ax.set_ylim(0, 12)
ax.set_xlabel(r"$\varepsilon_{ee}^{dV}$")
ax.set_ylabel(r"$\Delta\chi^{2}$")
fig.tight_layout(); fig.savefig(f"{BASE}/fig_conus_perfil_eps.png"); plt.close(fig)
print(f"  -> fig_conus_espectro.png, fig_conus_perfil_eps.png   (A_best={A}; 2a solucion en eps={ux[i2]:.3f}, dchi2={d2[i2]:.2f})")
