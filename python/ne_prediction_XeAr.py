#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prediccion SM del numero de eventos CEvNS en electrones de ionizacion extraidos,
Xe y Ar superpuestos, a 192 kg*dia (comparacion ideal, F de NEST). Marcadores
rellenos = dentro de la ROI de cada blanco (Xe N_e=4-7, Ar N_e=1-5); huecos = fuera.

Entrada: datos/espectro_Ne_ideal_{Xe,Ar}.dat (col 2 = R_bin con F de NEST)
Salida : datos/fig_ne_prediction_XeAr.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from estilo_tesis import aplicar, C_XE, C_AR, FIG15
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
EXPO = 192.0
ROI = {"Xe": (4, 7), "Ar": (1, 5)}
KMAX = {"Xe": 10, "Ar": 12}

fig, ax = plt.subplots(figsize=FIG15)
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    d = np.loadtxt(f"{BASE}/espectro_Ne_ideal_{tag}.dat", comments="#")
    k, N = d[:, 0].astype(int), d[:, 1] * EXPO
    sel = (k <= KMAX[tag]) & (N > 0)
    k, N = k[sel], N[sel]
    lo, hi = ROI[tag]
    en = (k >= lo) & (k <= hi)
    ax.plot(k, N, color=col, lw=1.6, alpha=0.7)
    ax.plot(k[en], N[en], "o", color=col, ms=8)
    ax.plot(k[~en], N[~en], "o", mfc="white", mec=col, mew=1.6, ms=8)

ax.set_yscale("log")
ax.set_xlim(0.5, 12.5)
ax.set_xticks(range(1, 13))
ax.set_xlabel(r"$N_e$  (electrones extraídos)")
ax.set_ylabel(r"eventos esperados en 192 kg$\cdot$día")
ax.legend(handles=[
    Line2D([], [], color=C_XE, marker="o", label="Xe"),
    Line2D([], [], color=C_AR, marker="o", label="Ar"),
    Line2D([], [], color="0.3", marker="o", ls="", label="dentro de la ROI"),
    Line2D([], [], color="0.3", marker="o", mfc="white", mew=1.6, ls="", label="fuera de la ROI")],
    loc="lower left")
fig.tight_layout()
out = f"{BASE}/fig_ne_prediction_XeAr.png"
fig.savefig(out); print(f"  -> {out}")
