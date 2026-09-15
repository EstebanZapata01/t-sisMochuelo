#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prediccion del modelo estandar del numero de eventos CEvNS en electrones de
ionizacion extraidos, Xe y Ar, a la exposicion de referencia 192 kg*dia
(comparacion ideal, F de NEST). Lee datos/espectro_Ne_ideal_{Xe,Ar}.dat
(col 2 = R_bin con la fluctuacion de NEST).

El ROI de cada blanco va sombreado:
  Xe : N_e = 4-5 ... 7  (arXiv:2411.18641; el piso N_e>=4 es fondo de SE)
  Ar : N_e = 1-5        (ref.[46] / SVII de arXiv:2411.18641)

Salida: datos/fig_ne_prediction_XeAr.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
EXPO = 192.0                       # kg*dia, tabla tab:ideal
ROI = {"Xe": (4, 7), "Ar": (1, 5)}
KMAX = {"Xe": 10, "Ar": 15}

fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), sharey=False)

for ax, (tag, col) in zip(axes, (("Xe", C_XE), ("Ar", C_AR))):
    d = np.loadtxt(f"{BASE}/espectro_Ne_ideal_{tag}.dat", comments="#")
    k, R = d[:, 0].astype(int), d[:, 1]
    N = R * EXPO                                   # eventos esperados
    lo, hi = ROI[tag]

    sel = k <= KMAX[tag]
    ax.bar(k[sel], np.maximum(N[sel], 1e-30), width=0.82, color=col,
           alpha=0.30, edgecolor=col, linewidth=1.0)
    ax.axvspan(lo - 0.5, hi + 0.5, color="0.45", alpha=0.14, lw=0,
               label=rf"ROI  $N_e = {lo}$–${hi}$")

    ax.set_yscale("log")
    ax.set_xlim(0.4, KMAX[tag] + 0.6)
    ax.set_xticks(range(1, KMAX[tag] + 1))
    ax.set_xlabel(r"$N_e$  (electrones extraídos)")
    ax.set_title(tag)
    ax.legend(loc="upper right")

axes[0].set_ylabel(r"eventos esperados en 192 kg$\cdot$día")
fig.suptitle(r"Predicción SM del espectro en $N_e$ (192 kg$\cdot$día)", y=1.02)
fig.tight_layout()
out = f"{BASE}/fig_ne_prediction_XeAr.png"
fig.savefig(out)
print(f"  -> {out}")
