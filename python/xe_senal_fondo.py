#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe (RED-100): senal CEvNS esperada y fondo medido (reactor apagado) por bin de N_e reconstruido, antes y
despues de los cortes de seleccion. Muestra que el fondo supera a la senal en ~3 ordenes de magnitud aun
despues de los cortes.

Entrada: datos/red100_senal_fondo_fig6.csv (python/extrae_pdf_red100.py)   Salida: datos/fig_xe_senal_fondo.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, ROJO, FIG15
aplicar(grande=True)

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
ne, sa, sd, fa, fd = np.loadtxt(f"{ROOT}/datos/red100_senal_fondo_fig6.csv", delimiter=",", comments="#", skiprows=4, unpack=True)

fig, ax = plt.subplots(figsize=FIG15)
for y, col, mk, hueco, lab in ((fa, ROJO, "s", True, "fondo, antes de cortes"), (fd, ROJO, "s", False, "fondo, después de cortes"),
                               (sa, C_XE, "o", True, "señal, antes de cortes"), (sd, C_XE, "o", False, "señal, después de cortes")):
    ax.plot(ne, y, color=col, alpha=0.3, lw=1.2)
    ax.plot(ne, y, mk, color=col, mfc="white" if hueco else col, mew=1.8, ms=9, label=lab)
ax.set_yscale("log")
ax.set_xticks(ne); ax.set_xticklabels([f"{int(n)}" for n in ne])
ax.set_xlim(3.5, 7.5); ax.set_ylim(1e-4, 1e4)
ax.set_xlabel(r"$N_e$ reconstruido")
ax.set_ylabel(r"cuentas / (kg$\cdot$día)")
ax.legend(loc="upper right", fontsize=10.5)
fig.tight_layout()
fig.savefig(f"{ROOT}/datos/fig_xe_senal_fondo.png"); print("  -> datos/fig_xe_senal_fondo.png")
