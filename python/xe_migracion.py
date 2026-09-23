#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe real: de que N_e VERDADERO (k) proviene cada bin de N_e RECONSTRUIDO (j=4..7),
en porcentaje del bin. Ilustra la migracion por resolucion en PE.

Entrada: datos/migracion_Xe.dat (red100PE.f90; R_k*P(k->j))
Salida : datos/fig_xe_migracion.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, GRIS, FIG15
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
j, k, c = np.loadtxt(f"{BASE}/migracion_Xe.dat", comments="#", unpack=True)
j, k = j.astype(int), k.astype(int)
BINS = [4, 5, 6, 7]
KS = [3, 4, 5, 6, 7]
COL = ["#c6d4df", "#8fa9bc", "#5c7f99", "#33546e", "#1f3547"]   # azules de menor a mayor k

tot = {b: c[j == b].sum() for b in BINS}
frac = {kk: np.array([100 * c[(j == b) & (k == kk)].sum() / tot[b] for b in BINS]) for kk in KS}
otros = np.array([100 - sum(frac[kk][i] for kk in KS) for i in range(len(BINS))])

fig, ax = plt.subplots(figsize=FIG15)
base = np.zeros(len(BINS))
for kk, col in zip(KS, COL):
    ax.bar(BINS, frac[kk], bottom=base, color=col, width=0.62, label=rf"$k={kk}$")
    for x, h, b0 in zip(BINS, frac[kk], base):
        if h >= 6:
            ax.text(x, b0 + h / 2, f"{h:.0f}%", ha="center", va="center", fontsize=11,
                    color="white" if col in COL[2:] else "0.15")
    base += frac[kk]
ax.bar(BINS, otros, bottom=base, color=GRIS, width=0.62, alpha=0.6, label="otros")
ax.set_xticks(BINS)
ax.set_xlabel(r"$N_e$ reconstruido  ($j$)")
ax.set_ylabel("procedencia [% del bin]")
ax.set_ylim(0, 100)
ax.legend(title=r"$N_e$ verdadero", loc="center left", bbox_to_anchor=(1.01, 0.5))
fig.tight_layout()
out = f"{BASE}/fig_xe_migracion.png"
fig.savefig(out); print(f"  -> {out}")
