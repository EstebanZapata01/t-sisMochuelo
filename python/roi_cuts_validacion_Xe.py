#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe, ROI (N_e = 4-7), todo en N_e RECONSTRUIDO (PE/27): senal de RED-100
antes/despues de cortes (digitalizada del paper, arXiv:2411.18641) frente a
la simulacion pasada por la misma respuesta en PE y la misma ventana.

Entrada: datos/validacion_fig3_fig6_Xe.dat (red100PE.f90; la simulacion reconstruida
         sale de prob_migracion en mod_detector.f90)
Salida : datos/fig_roi_cuts_Xe.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, NEGRO
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0, "lines.linewidth": 1.3,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"

d = np.loadtxt(f"{BASE}/validacion_fig3_fig6_Xe.dat", comments="#")
# Ne sim_verdadero sim_reconstruido sim_rec_x_effROI fig3_extraidos fig6_antes fig6_despues
Ne, sim, ajuste = d[:, 0], d[:, 2], d[:, 3]
red_antes, red_despues = d[:, 5], d[:, 6]

fig, ax = plt.subplots(figsize=(4.91, 3.6))
series = [                     # (datos, color, marcador, etiqueta, relleno)
    (red_antes,   NEGRO, "o", "publicado, antes de cortes"),
    (red_despues, NEGRO, "^", "publicado, después de cortes"),
    (sim,         C_XE,  "o", "simulación, antes de cortes"),
    (ajuste,      C_XE,  "^", "simulación, después de cortes"),
]
for y, col, mk, lab in series:
    ax.plot(Ne, y, color=col, alpha=0.25, lw=1.2, zorder=2)
    ax.plot(Ne, y, mk, color=col, ms=7, label=lab, zorder=3)
ax.set_yscale("log")
ax.set_xticks(Ne); ax.set_xticklabels([f"{int(n)}" for n in Ne])
ax.set_xlim(3.5, 7.5); ax.set_ylim(1.2e-4, 3e-2)
ax.set_xlabel(r"$N_e$ reconstruido (PE$/27$)")
ax.set_ylabel(r"eventos / (kg$\cdot$día)")
ax.legend(loc="upper right", fontsize=8.0)
fig.tight_layout()
out = f"{BASE}/fig_roi_cuts_Xe.png"
fig.savefig(out, dpi=300)
dest = f"{DEST}/fig_roi_cuts_Xe.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print("Ne  sim/RED-100(antes)  ajuste/RED-100(despues)")
for i in range(len(Ne)):
    print(f"{int(Ne[i])}   {sim[i]/red_antes[i]:.2f}   {ajuste[i]/red_despues[i]:.2f}")
print(f"\n  -> {dest}  (md5 {md5(out)})")
