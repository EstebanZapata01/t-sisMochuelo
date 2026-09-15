#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Efecto MEDIDO de la aproximacion F^2 = 1, desde el motor real.

Corre red100_nest (Xe y Ar) dos veces cada uno --- sin y con el factor de
forma de Helm activado (variable de entorno USE_HELM, leida por
xsections_nest.f90) --- y grafica el cociente

    R(T_nr) = [dR/dT_nr]_Helm  /  [dR/dT_nr]_{F^2=1}

que es exactamente 1 - F^2_Helm(q^2(T_nr)) reproducido por la simulacion.
Sombra: el tramo de T_nr dentro del ROI de ionizacion de cada blanco.

Restaura los .dat de la linea base (F^2 = 1) al terminar.

Entrada/salida por el motor:  datos/espectro_continuo{Xe,Ar}.dat
Salida figura:  datos/fig_helm_FF.png
"""
import os
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR
aplicar()

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DIR = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe",
       "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}
T_THR = {"Xe": 0.64, "Ar": 0.26}          # keV_nr, borde inferior del ROI
BK = "/tmp/helm_ff_baseline"


def run_red100(tag, use_helm):
    e = dict(os.environ)
    if use_helm:
        e["USE_HELM"] = "1"
    else:
        e.pop("USE_HELM", None)
    subprocess.run(["./red100_nest"], cwd=DIR[tag], env=e, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    d = np.loadtxt(f"{BASE}/espectro_continuo{tag}.dat", comments="#")
    return d[:, 0], d[:, 3]                # T_nr, dR/dT_nr (Comb, hibrido)


os.makedirs(BK, exist_ok=True)
for tag in ("Xe", "Ar"):
    shutil.copy(f"{BASE}/espectro_continuo{tag}.dat", f"{BK}/{tag}.dat")

fig, ax = plt.subplots(figsize=(6.8, 4.4))
try:
    for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
        T0, r0 = run_red100(tag, use_helm=False)      # F^2 = 1
        T1, r1 = run_red100(tag, use_helm=True)       # Helm
        m = (r0 > 0) & (T0 <= 4.0)
        ratio = np.where(m, r1 / np.where(r0 == 0, np.nan, r0), np.nan)
        ax.plot(T0[m], ratio[m], color=col, lw=1.8, label=tag)
        ax.axvline(T_THR[tag], color=col, ls=":", lw=1.1)
        ax.text(T_THR[tag] + 0.04, 0.30, rf"$T_{{\rm thr}}^{{\rm {tag}}}$",
                rotation=90, va="bottom", ha="left", fontsize=8, color=col,
                transform=ax.get_xaxis_transform())
finally:
    for tag in ("Xe", "Ar"):
        shutil.copy(f"{BK}/{tag}.dat", f"{BASE}/espectro_continuo{tag}.dat")
        run_red100(tag, use_helm=False)               # deja el .dat en F^2=1

ax.axhline(1.0, color="0.6", lw=0.9, ls="--")
ax.set_xlim(0.2, 4.0)
ax.set_ylim(0.90, 1.005)
ax.set_xlabel(r"Energía de retroceso nuclear  $T_{\rm nr}$  [keV]")
ax.set_ylabel(r"$(dR/dT_{\rm nr})_{\rm Helm}\,/\,(dR/dT_{\rm nr})_{F^2=1}$")
ax.set_title("Efecto medido de la aproximación $F^2=1$ (motor: red100_nest, USE_HELM)")
ax.legend(loc="center right", title="CE$\\nu$NS SM")

fig.tight_layout()
out = f"{BASE}/fig_helm_FF.png"
fig.savefig(out)
print(f"  -> {out}")
