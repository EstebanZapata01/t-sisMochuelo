#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Argon: senal CEvNS S_k y fondo de 39Ar B_k (atmosferico y depletado) por bin de
N_e = 1-5 (la ROI de Ar), a 62 kg*dia. Muestra que el fondo atmosferico
domina desde N_e ~ 4 y que el depletado es despreciable.

Entrada: datos/bkg_bins_Ar.dat (chi2_bkg_nest.f90; ev/(kg dia))
Salida : datos/fig_ar_senal_fondo.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_AR, ROJO, MORADO, FIG15
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
EXPO = 62.0                                    # kg*dia (masa de la ref. de argon, 1 dia)
k, S, Bu, Ba = np.loadtxt(f"{BASE}/bkg_bins_Ar.dat", comments="#", unpack=True)

fig, ax = plt.subplots(figsize=FIG15)
ax.plot(k, S * EXPO, "o-", color=C_AR, label="señal CE$\\nu$NS")
ax.plot(k, Ba * EXPO, "s-", color=ROJO, label=r"fondo $^{39}$Ar atmosférico")
ax.plot(k, Bu * EXPO, "^-", color=MORADO, label=r"fondo $^{39}$Ar depletado")
ax.set_yscale("log")
ax.set_xticks([1, 2, 3, 4, 5]); ax.set_xlim(0.6, 5.4)
ax.set_xlabel(r"$N_e$  (electrones extraídos)")
ax.set_ylabel(r"eventos en 62 kg$\cdot$día")
ax.legend(loc="lower left")
fig.tight_layout()
out = f"{BASE}/fig_ar_senal_fondo.png"
fig.savefig(out); print(f"  -> {out}")
