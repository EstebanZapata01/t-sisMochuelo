#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Espectro de antineutrinos de reactor usado (hibrido: Kopeikin por debajo de 2 MeV,
Huber-Mueller desde 2 MeV; normalizado a 1.4e13 cm^-2 s^-1), con la energia
minima E_nu^min = sqrt(M T/2) para T_nr = 0.2 keV (suelo de NEST) en Xe y Ar.

Entrada: datos_tutor/insumos/flujo_nu_hibrido.csv (exporta_datos_tutor.py, con flux.f90)
Salida : datos/fig_flujo_antineutrinos.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, NEGRO, FIG15
aplicar(grande=True)

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
E, F = np.loadtxt(f"{ROOT}/datos_tutor/insumos/flujo_nu_hibrido.csv", delimiter=",", skiprows=1, unpack=True)
AMU = 931.49410242                                        # MeV/u
E_MIN = {"Xe": np.sqrt(131.293 * AMU * 2e-4 / 2), "Ar": np.sqrt(39.948 * AMU * 2e-4 / 2)}

fig, ax = plt.subplots(figsize=FIG15)
sel = F > 0
ax.plot(E[sel], F[sel], color=NEGRO)
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    ax.axvline(E_MIN[tag], color=col, ls="--", lw=1.4, label=rf"$E_\nu^{{\min}}$ ({tag}, 0,2 keV)")
ax.set_yscale("log")
ax.set_xlim(0, 10); ax.set_ylim(1e6, 2e13)
ax.set_xlabel(r"$E_\nu$  [MeV]")
ax.set_ylabel(r"$d\Phi/dE_\nu$  [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]")
ax.legend(loc="upper right")
fig.tight_layout()
out = f"{ROOT}/datos/fig_flujo_antineutrinos.png"
fig.savefig(out); print(f"  -> {out}")
