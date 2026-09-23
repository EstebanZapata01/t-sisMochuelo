#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Espectro en electrones de ionizacion (creados y extraidos) de Xe: simulacion
frente al espectro publicado por RED-100 (digitalizado).

Entradas: datos/ionization_electrones.dat (red100_nest.f90), datos/red100_fig3_digitalizado.csv
Salida  : datos/fig_ionizationXe.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, NEGRO, FIG15
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
sim = np.loadtxt(f"{BASE}/ionization_electrones.dat", comments="#")            # Ne creados extraidos
pub = np.loadtxt(f"{BASE}/red100_fig3_digitalizado.csv", delimiter=",", skiprows=1)
# Ne = 0 (extraidos publicados) no se usa: el punto digitalizado no conserva el numero de eventos.
sim, pub = sim[(sim[:, 0] >= 1) & (sim[:, 0] <= 10)], pub[pub[:, 0] >= 1]

fig, ax = plt.subplots(figsize=FIG15)
for col, mk, nombre in ((1, "o", "creados"), (2, "s", "extraídos")):
    ax.plot(sim[:, 0], sim[:, col], color=C_XE, lw=1.1, alpha=0.5)
    ax.plot(sim[:, 0], sim[:, col], mk, color=C_XE, ms=7, label=f"{nombre}, simulación")
    ax.plot(pub[:, 0], pub[:, col], mk, color=NEGRO, mfc="white", mew=1.4, ms=7, label=f"{nombre}, publicado")
ax.set_yscale("log")
ax.set_xlim(0.5, 10.5); ax.set_xticks(range(1, 11))
ax.set_ylim(1e-8, 1e2)
ax.set_xlabel(r"$N_e$  (electrones de ionización)")
ax.set_ylabel(r"eventos / (kg $\cdot$ día)")
ax.legend(loc="upper right", ncol=1)
fig.tight_layout()
out = f"{BASE}/fig_ionizationXe.png"
fig.savefig(out); print(f"  -> {out}")
