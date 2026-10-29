#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Espectro en electrones de ionizacion (creados y extraidos) de Xe: simulacion
frente al espectro publicado por RED-100 (digitalizado).

Entradas: datos/ionization_electrones.dat (red100_nest.f90), datos/red100_fig3_digitalizado.csv
Salida  : datos/fig_ionizationXe.png (+ copia en figuras_metodologia)
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
sim = np.loadtxt(f"{BASE}/ionization_electrones.dat", comments="#")            # Ne creados extraidos
pub = np.loadtxt(f"{BASE}/red100_fig3_digitalizado.csv", delimiter=",", skiprows=1)
# Ne = 0 (extraidos publicados) no se usa: el punto digitalizado no conserva el numero de eventos.
sim, pub = sim[(sim[:, 0] >= 1) & (sim[:, 0] <= 10)], pub[pub[:, 0] >= 1]

fig, ax = plt.subplots(figsize=(4.91, 3.6))
for col, mk, nombre in ((1, "o", "creados"), (2, "s", "extraídos")):
    ax.plot(sim[:, 0], sim[:, col], color=C_XE, lw=1.1, alpha=0.5)
    ax.plot(sim[:, 0], sim[:, col], mk, color=C_XE, ms=6, label=f"{nombre}, simulación")
    ax.plot(pub[:, 0], pub[:, col], mk, color=NEGRO, mfc="white", mew=1.3, ms=6, label=f"{nombre}, publicado")
ax.set_yscale("log")
ax.set_xlim(0.5, 10.5); ax.set_xticks(range(1, 11))
ax.set_ylim(1e-8, 1e2)
ax.set_xlabel(r"$N_e$  (electrones de ionización)")
ax.set_ylabel(r"eventos / (kg $\cdot$ día)")
ax.legend(loc="upper right", ncol=1)
fig.tight_layout()
out = f"{BASE}/fig_ionizationXe.png"
fig.savefig(out, dpi=300)
dest = f"{DEST}/fig_ionizationXe.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"  -> {dest}  (md5 {md5(out)})")
