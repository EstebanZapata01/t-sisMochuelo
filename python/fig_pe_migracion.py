#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe: espectro simulado en energia corregida [PE] y migracion hacia la ROI
reconstruida. Fusiona fig_espectro_PE_Xe.py (arriba) y xe_migracion.py
(abajo) en una sola figura, 0,9 del ancho de pagina.

Arriba: espectro completo (0-220 PE), total y plantillas k=1..7 (campana
N(27k, sqrt(k) 7,6)), ROI 110-189 PE sombreada, bordes de ventana j=4..7
marcados.
Abajo: procedencia por bin reconstruido j=4..7 (igual que xe_migracion.py):
barras cuyo ancho es la ventana real en PE (recortada a la ROI) y cuya
altura apilada por k es el % del area de esa ventana que aporta cada k.

Entrada: datos/ionization_spectra_detallado.dat, datos/migracion_Xe.dat
         (red100PE.f90)
Salida : datos/fig_pe_migracion.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from estilo_tesis import aplicar, GRIS, NEGRO
from leer_fortran import leer_dat
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8, "lines.linewidth": 1.15,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
ROI_LO, ROI_HI = 110.0, 189.0
EDGES = [27.0 * (j - 0.5) for j in (4, 5, 6, 7, 8)]   # 94.5..202.5
KS = [3, 4, 5, 6, 7]
COLK = ["#c6d4df", "#8fa9bc", "#5c7f99", "#33546e", "#1f3547"]
AZULES7 = ["#dbe4ea", "#c6d4df", "#a3b9ca", "#7d9bb2", "#5c7f99", "#33546e", "#1f3547"]

d, _ = leer_dat(f"{BASE}/ionization_spectra_detallado.dat")   # PE Total 1SE..7SE
pe, tot = d[:, 0], d[:, 1]

jj, kk, cc = np.loadtxt(f"{BASE}/migracion_Xe.dat", comments="#", unpack=True)
jj, kk = jj.astype(int), kk.astype(int)
BINS = [4, 5, 6, 7]
tot_j = {b: cc[jj == b].sum() for b in BINS}
frac = {k: np.array([100 * cc[(jj == b) & (kk == k)].sum() / tot_j[b] for b in BINS]) for k in KS}
otros = np.array([100 - sum(frac[k][i] for k in KS) for i in range(len(BINS))])

fig, (a1, a2) = plt.subplots(2, 1, figsize=(5.53, 5.0),
                              gridspec_kw=dict(height_ratios=[1.6, 1.0], hspace=0.30))

# ------------------------------------------------------------- arriba: espectro completo
a1.axvspan(ROI_LO, ROI_HI, color=GRIS, alpha=0.12, lw=0)
a1.fill_between(pe, 1e-30, tot, step="mid", color=GRIS, alpha=0.25, lw=0, label="total")
for k in range(1, 8):
    a1.plot(pe, d[:, 1 + k], color=AZULES7[k - 1], lw=1.3, label=rf"$k={k}$")
for edge in EDGES:
    a1.axvline(edge, color="0.5", lw=0.7, ls=(0, (4, 3)))
trans = mtransforms.blended_transform_factory(a1.transData, a1.transAxes)
for j, elo, ehi in zip((4, 5, 6, 7), EDGES[:-1], EDGES[1:]):
    a1.text(0.5 * (elo + ehi), 0.58, f"$j={j}$", ha="center", va="top",
            fontsize=7.5, color=NEGRO, transform=trans)
a1.set_yscale("log")
a1.set_xlim(0, 225); a1.set_ylim(1e-5, 8)
a1.set_xlabel("energía corregida $S$ [PE]")
a1.set_ylabel("eventos /\n(5 PE $\\cdot$ kg $\\cdot$ día)")
a1.legend(loc="upper right", fontsize=6.6, ncol=4, columnspacing=0.8, handlelength=1.3)

# ------------------------------------------------------------- abajo: procedencia por bin
base = np.zeros(len(BINS))
win_lo = [max(e, ROI_LO) for e in EDGES[:-1]]
win_hi = [min(e, ROI_HI) for e in EDGES[1:]]
for k, col in zip(KS, COLK):
    for i, b in enumerate(BINS):
        h = frac[k][i]
        a2.bar(win_lo[i], h, bottom=base[i], width=win_hi[i] - win_lo[i],
               align="edge", color=col, edgecolor="white", lw=0.5)
        if h >= 10:
            a2.text(0.5 * (win_lo[i] + win_hi[i]), base[i] + h / 2, f"$k={k}$ · {h:.0f}%",
                    ha="center", va="center", fontsize=6.8,
                    color="white" if col in COLK[2:] else "0.15")
    base += frac[k]
for i, b in enumerate(BINS):
    a2.bar(win_lo[i], otros[i], bottom=base[i], width=win_hi[i] - win_lo[i],
           align="edge", color=GRIS, alpha=0.6, edgecolor="white", lw=0.5)
a2.set_ylim(0, 105)
a2.set_xlim(108, 191)
a2.set_xlabel("energía corregida $S$ [PE]")
a2.set_ylabel("% de la ventana")

fig.tight_layout()
out = f"{BASE}/fig_pe_migracion.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = f"{DEST}/fig_pe_migracion.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
