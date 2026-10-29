#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prediccion SM del numero de eventos CEvNS en electrones de ionizacion extraidos,
Xe y Ar superpuestos, a 192 kg*dia (comparacion ideal, F de NEST). Marcadores
rellenos = dentro de la ROI de cada blanco (Xe N_e=4-7, Ar N_e=1-5); huecos = fuera.

Entrada: datos/espectro_Ne_ideal_{Xe,Ar}.dat (col 2 = R_bin con F de NEST)
Salida : datos/fig_ne_prediction_XeAr.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from estilo_tesis import aplicar, C_XE, C_AR
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0, "lines.linewidth": 1.3,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
EXPO = 192.0
ROI = {"Xe": (4, 7), "Ar": (1, 5)}
KMAX = {"Xe": 10, "Ar": 12}

fig, ax = plt.subplots(figsize=(4.91, 3.6))
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    d = np.loadtxt(f"{BASE}/espectro_Ne_ideal_{tag}.dat", comments="#")
    k, N = d[:, 0].astype(int), d[:, 1] * EXPO
    sel = (k <= KMAX[tag]) & (N > 0)
    k, N = k[sel], N[sel]
    lo, hi = ROI[tag]
    en = (k >= lo) & (k <= hi)
    ax.plot(k, N, color=col, lw=1.4, alpha=0.7)
    ax.plot(k[en], N[en], "o", color=col, ms=6.5)
    ax.plot(k[~en], N[~en], "o", mfc="white", mec=col, mew=1.4, ms=6.5)

ax.set_yscale("log")
ax.set_xlim(0.5, 12.5)
ax.set_xticks(range(1, 13))
ax.set_xlabel(r"$N_e$  (electrones extraídos)")
ax.set_ylabel(r"eventos esperados en 192 kg$\cdot$día")
ax.legend(handles=[
    Line2D([], [], color=C_XE, marker="o", label="Xe"),
    Line2D([], [], color=C_AR, marker="o", label="Ar"),
    Line2D([], [], color="0.3", marker="o", ls="", label="dentro de la ROI"),
    Line2D([], [], color="0.3", marker="o", mfc="white", mew=1.4, ls="", label="fuera de la ROI")],
    loc="lower left")
fig.tight_layout()
out = f"{BASE}/fig_ne_prediction_XeAr.png"
fig.savefig(out, dpi=300)
dest = f"{DEST}/fig_ne_prediction_XeAr.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"  -> {dest}  (md5 {md5(out)})")
