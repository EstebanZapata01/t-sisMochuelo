#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe: espectro simulado en energia corregida [PE]: total y plantillas de 1 a 6
electrones (campana N(27 k, sqrt(k) 7.6)), con la ROI 110-189 PE sombreada.

Entrada: datos/ionization_spectra_detallado.dat (red100PE.f90; ev/(5 PE kg dia))
Salida : datos/fig_espectro_PE_Xe.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, GRIS, NEGRO, FIG15
from leer_fortran import leer_dat
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
d, _ = leer_dat(f"{BASE}/ionization_spectra_detallado.dat")     # PE Total 1SE..7SE
pe, tot = d[:, 0], d[:, 1]
AZULES = ["#c6d4df", "#a3b9ca", "#7d9bb2", "#5c7f99", "#33546e", "#1f3547"]

fig, ax = plt.subplots(figsize=FIG15)
ax.axvspan(110, 189, color=GRIS, alpha=0.12, lw=0)
ax.fill_between(pe, 1e-30, tot, step="mid", color=GRIS, alpha=0.25, lw=0, label="total")
for k in range(1, 7):
    ax.plot(pe, d[:, 1 + k], color=AZULES[k - 1], lw=1.9, label=rf"$k={k}$")
ax.set_yscale("log")
ax.set_xlim(0, 225); ax.set_ylim(1e-5, 5)
ax.set_xlabel("energía corregida [PE]")
ax.set_ylabel(r"eventos / (5 PE $\cdot$ kg $\cdot$ día)")
ax.legend(loc="upper right", ncol=2)
fig.tight_layout()
out = f"{BASE}/fig_espectro_PE_Xe.png"
fig.savefig(out); print(f"  -> {out}")
