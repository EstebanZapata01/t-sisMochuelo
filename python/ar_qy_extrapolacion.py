#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Argon: Q_y(T_nr) de LArNEST reescalado para pasar por ReD, con el rango medido
(energias medias de los puntos, 2.4-7.6 keV) y la ROI CEvNS (0.1-1 keV)
sombreada: lo que queda por debajo de 2.4 keV es extrapolacion.
La tabla llega a 6 keV, asi que el rango medido se dibuja hasta el borde.

Entrada: datos/nest_Ar_218V_dense.txt   Salida: datos/fig_ar_qy_extrapolacion.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_AR, GRIS, FIG15
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
T, Qy = np.loadtxt(f"{BASE}/nest_Ar_218V_dense.txt", comments="#", usecols=(0, 1), unpack=True)
T_MED_MIN = 2.4    # keV: energia media del primer punto medido por ReD

fig, ax = plt.subplots(figsize=FIG15)
ax.axvspan(0.1, 1.0, color=GRIS, alpha=0.20, lw=0)
ax.axvspan(T_MED_MIN, 6.0, color=C_AR, alpha=0.12, lw=0)
ext = T < T_MED_MIN
ax.plot(T[ext], Qy[ext], color=C_AR, ls="--", label="extrapolación")
ax.plot(T[~ext], Qy[~ext], color=C_AR, label="rango medido (2,4–7,6 keV)")
ax.text(0.55, 0.05, "ROI\nCE$\\nu$NS", transform=ax.get_xaxis_transform(), ha="center",
        va="bottom", fontsize=11, color="0.35")
ax.set_xlim(0, 6)
ax.set_ylim(0, None)
ax.set_xlabel(r"$T_{\rm nr}$  [keV]")
ax.set_ylabel(r"$Q_y(T_{\rm nr})$  [e$^-$/keV]")
ax.legend(loc="lower right")
fig.tight_layout()
out = f"{BASE}/fig_ar_qy_extrapolacion.png"
fig.savefig(out); print(f"  -> {out}")
