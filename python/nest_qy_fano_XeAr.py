#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rendimiento de carga Q_y(T) y factor de Fano F(T)=Var(N_e)/<N_e> de NEST (Xe)
y LArNEST (Ar), superpuestos, hasta T~2 keV (la ROI CEvNS de ambos blancos).

Entradas: datos/nest_218V_dense.txt, datos/nest_Ar_218V_dense.txt
Salidas : datos/fig_nest_qy.png, datos/fig_nest_fano.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, FIG15
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
T_MAX = 2.0   # keV
dXe = np.loadtxt(f"{BASE}/nest_218V_dense.txt", comments="#")
dAr = np.loadtxt(f"{BASE}/nest_Ar_218V_dense.txt", comments="#")

for col, nombre, ylab, ylim, loc in (
        (1, "fig_nest_qy", r"$Q_y(T_{\rm nr})$  [e$^-$/keV]", (0, None), "lower right"),
        (2, "fig_nest_fano", r"$F(T_{\rm nr}) = \mathrm{Var}(N_e)/\langle N_e\rangle$", (0, 0.7), "center right")):
    fig, ax = plt.subplots(figsize=FIG15)
    for d, c, tag in ((dXe, C_XE, "Xe"), (dAr, C_AR, "Ar")):
        sel = d[:, 0] <= T_MAX
        ax.plot(d[sel, 0], d[sel, col], color=c, label=tag)
    ax.set_xlim(0, T_MAX); ax.set_ylim(*ylim)
    ax.set_xlabel(r"$T_{\rm nr}$  [keV]"); ax.set_ylabel(ylab)
    ax.legend(loc=loc)
    fig.tight_layout(); fig.savefig(f"{BASE}/{nombre}.png"); plt.close(fig)
    print(f"  -> {BASE}/{nombre}.png")
