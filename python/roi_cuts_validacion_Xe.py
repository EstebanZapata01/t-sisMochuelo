#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe, ROI (N_e = 4-7), todo en N_e RECONSTRUIDO (PE/27): senal de RED-100
antes/despues de cortes (digitalizada del paper, arXiv:2411.18641) frente a
la simulacion pasada por la misma respuesta en PE y la misma ventana.

Entrada: datos/validacion_fig3_fig6_Xe.dat (red100PE.f90; la simulacion reconstruida
         sale de prob_migracion en mod_detector.f90)
Salida : datos/fig_roi_cuts_Xe.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, NEGRO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

d = np.loadtxt(f"{BASE}/validacion_fig3_fig6_Xe.dat", comments="#")
# Ne sim_verdadero sim_reconstruido sim_rec_x_effROI fig3_extraidos fig6_antes fig6_despues
Ne, sim, ajuste = d[:, 0], d[:, 2], d[:, 3]
red_antes, red_despues = d[:, 5], d[:, 6]

fig, ax = plt.subplots(figsize=(7.0, 5.2))
series = [
    (red_antes,   NEGRO, "^", "RED-100, antes de cortes"),
    (red_despues, NEGRO, "s", "RED-100, después de cortes"),
    (sim,         C_XE,  "o", "Simulación"),
    (ajuste,      C_XE,  "s", r"Simulación $\times\,\varepsilon_{\rm ROI}$ (ajuste)"),
]
for y, col, mk, lab in series:
    ax.plot(Ne, y, color=col, alpha=0.3, lw=1.2, zorder=2)
    ax.scatter(Ne, y, color=col, marker=mk, s=42, label=lab, zorder=3)

ax.set_yscale("log")
ax.set_xticks(Ne); ax.set_xticklabels([f"{int(n)}" for n in Ne])
ax.set_xlim(3.5, 7.5)
ax.set_xlabel(r"$N_e$ reconstruido (PE$/27$)")
ax.set_ylabel(r"eventos / (kg$\cdot$día)")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
out = f"{BASE}/fig_roi_cuts_Xe.png"
fig.savefig(out)

print("Ne  sim/RED-100(antes)  ajuste/RED-100(despues)")
for i in range(len(Ne)):
    print(f"{int(Ne[i])}   {sim[i]/red_antes[i]:.2f}   {ajuste[i]/red_despues[i]:.2f}")
print(f"\n  {out}")
