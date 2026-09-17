#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cinematica exacta CEvNS: T_max(E_nu) = 2 E_nu^2 / (M + 2 E_nu), sin la
aproximacion E_nu << M. Xe (mezcla isotopica natural, <A>=131.293, la misma
que usa el resto del pipeline) vs Ar (A=40) -- la asimetria cinematica de la
que depende toda la comparacion Xe/Ar: a igual E_nu, el argon retrocede mas.

Masas: M = A * amu (mismos A_Ge, amu de constants.f90 en
N_EventosCEvNS_NSIXe/ y N_EventosCEvNS_NSIAr/; no se rederivan, se citan).

Salida: datos/fig_tmax_kinematica.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

AMU = 931.49410242          # MeV/u, constants.f90 (ambas carpetas)
A_XE, A_AR = 131.293, 39.948  # constants.f90: A_Ge (Xe natural, Ar natural)
M_XE, M_AR = A_XE * AMU, A_AR * AMU  # MeV


def Tmax_keV(E_nu_MeV, M_MeV):
    return 2.0 * E_nu_MeV**2 / (M_MeV + 2.0 * E_nu_MeV) * 1000.0


E_nu = np.linspace(1e-6, 10.0, 500)

fig, ax = plt.subplots(figsize=(6.8, 4.6))
ax.plot(E_nu, Tmax_keV(E_nu, M_XE), color=C_XE, label="Xe")
ax.plot(E_nu, Tmax_keV(E_nu, M_AR), color=C_AR, label="Ar")

ax.set_xlim(0, 10)
ax.set_ylim(0, None)
ax.set_xlabel(r"Energía del antineutrino, $E_\nu$  [MeV]")
ax.set_ylabel(r"$T_{\max}(E_\nu)$  [keV]")
ax.legend(loc="upper left")
fig.tight_layout()
out = f"{BASE}/fig_tmax_kinematica.png"
fig.savefig(out)

print("=" * 60)
print(f" M_Xe = {M_XE:.2f} MeV   M_Ar = {M_AR:.2f} MeV")
for E in (3.0, 5.0, 8.0):
    print(f"  E_nu={E:.0f} MeV: T_max Xe={Tmax_keV(E, M_XE):.4f} keV   "
          f"T_max Ar={Tmax_keV(E, M_AR):.4f} keV   "
          f"razon Ar/Xe={Tmax_keV(E, M_AR)/Tmax_keV(E, M_XE):.3f}")
print(f"\n  {out}")
