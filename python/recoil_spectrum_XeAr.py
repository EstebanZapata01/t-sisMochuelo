#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Espectro de retroceso nuclear CEvNS del modelo estandar, dR/dT_nr, para Xe y
Ar (flujo hibrido Kopeikin U Huber-Mueller). Es el punto de partida del
pipeline: lo que se integra contra la seccion eficaz antes de NEST.

Lee datos/espectro_continuo{Xe,Ar}.dat  (col 4 = Tasa_Comb, hibrido).
Marca T_thr de cada blanco (energia a la que <N_e> alcanza el borde del ROI,
tabla_umbral.tex).

Salida: datos/fig_recoil_XeAr.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
T_THR = {"Xe": 0.64, "Ar": 0.26}          # keV_nr, borde inferior del ROI

fig, ax = plt.subplots(figsize=(6.8, 4.4))

for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    d = np.loadtxt(f"{BASE}/espectro_continuo{tag}.dat", comments="#")
    T, dRdT = d[:, 0], d[:, 3]
    m = dRdT > 1e-3
    ax.plot(T[m], dRdT[m], color=col, lw=1.8, label=tag)
    ax.axvline(T_THR[tag], color=col, ls=":", lw=1.1)
    ax.text(T_THR[tag], 0.03, rf" $T_{{\rm thr}}^{{\rm {tag}}}={T_THR[tag]}$",
            rotation=90, va="bottom", ha="left", fontsize=8, color=col,
            transform=ax.get_xaxis_transform())

ax.set_yscale("log")
ax.set_xlim(0.2, 4.0)
ax.set_ylim(1e-3, 4e2)
ax.set_xlabel(r"Energía de retroceso nuclear  $T_{\rm nr}$  [keV]")
ax.set_ylabel(r"$dR/dT_{\rm nr}$  [ev$\cdot$kg$^{-1}\cdot$día$^{-1}\cdot$keV$^{-1}$]")
ax.legend(loc="upper right", title="CE$\\nu$NS SM")
ax.set_title(r"Espectro de retroceso nuclear CE$\nu$NS (SM)")

fig.tight_layout()
out = f"{BASE}/fig_recoil_XeAr.png"
fig.savefig(out)
print(f"  -> {out}")
