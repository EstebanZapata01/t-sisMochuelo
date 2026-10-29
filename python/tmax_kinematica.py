#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cinematica exacta CEvNS: T_max(E_nu) = 2 E_nu^2 / (M + 2 E_nu), sin la
aproximacion E_nu << M (esa aproximacion, 2E_nu^2/M, se dibuja punteada para
mostrar donde se separa de la formula exacta). Xe (mezcla isotopica natural,
<A>=131.293, la misma que usa el resto del pipeline) vs Ar (A=40) -- la
asimetria cinematica de la que depende toda la comparacion Xe/Ar: a igual
E_nu, el argon retrocede mas.

Umbral del analisis de RED-100 en Xe: N_e=4 equivale a T~2 keV de retroceso
nuclear (Sec. de electrones medios / fig_electrones_medios). Se marca esa
linea horizontal y se resalta donde el T_max exacto del Xe la cruza: dentro
del espectro del reactor (que se apaga fuerte por encima de 8 MeV), el Xe
casi nunca llega a esa T con un solo evento -- sus eventos en la ROI real
entran solo por la cola dura del flujo y por fluctuaciones NEST.

Masas: M = A * amu (mismos A_Ge, amu de constants.f90 en
N_EventosCEvNS_NSIXe/ y N_EventosCEvNS_NSIAr/; no se rederivan, se citan).

Salida: datos/fig_tmax_kinematica.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, NEGRO, GRIS, coma, formatear_coma
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0, "lines.linewidth": 1.3,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"

AMU = 931.49410242          # MeV/u, constants.f90 (ambas carpetas)
A_XE, A_AR = 131.293, 39.948  # constants.f90: A_Ge (Xe natural, Ar natural)
M_XE, M_AR = A_XE * AMU, A_AR * AMU  # MeV
T_UMBRAL = 2.0   # keV, ~4 electrones en Xe (fig_electrones_medios)
E_COLA = 8.0     # MeV, borde donde el espectro del reactor se apaga


def Tmax_keV(E_nu_MeV, M_MeV):
    return 2.0 * E_nu_MeV**2 / (M_MeV + 2.0 * E_nu_MeV) * 1000.0


def Tmax_aprox_keV(E_nu_MeV, M_MeV):
    return 2.0 * E_nu_MeV**2 / M_MeV * 1000.0


E_nu = np.linspace(1e-6, 12.0, 600)

fig, ax = plt.subplots(figsize=(4.60, 3.6))
ax.axvspan(E_COLA, 12.0, color=GRIS, alpha=0.15, zorder=0)
ax.text(E_COLA + 0.15, 0.94, "cola del\nespectro\ndel reactor",
        transform=ax.get_xaxis_transform(), ha="left", va="top", fontsize=8.0, color="0.35")

ax.plot(E_nu, Tmax_keV(E_nu, M_XE), color=C_XE, lw=1.8, label="Xe (exacto)")
ax.plot(E_nu, Tmax_keV(E_nu, M_AR), color=C_AR, lw=1.8, label="Ar (exacto)")
ax.plot(E_nu, Tmax_aprox_keV(E_nu, M_XE), color=C_XE, lw=1.2, ls=":", label=r"Xe ($2E_\nu^2/M$)")
ax.plot(E_nu, Tmax_aprox_keV(E_nu, M_AR), color=C_AR, lw=1.2, ls=":", label=r"Ar ($2E_\nu^2/M$)")

ax.axhline(T_UMBRAL, color=NEGRO, lw=1.1, ls="--")
ax.text(0.15, T_UMBRAL + 0.2, "umbral RED-100, Xe\n($\\approx4\\,e^-$)",
        fontsize=8.0, color=NEGRO, va="bottom")

E_cruce = brentq(lambda E: Tmax_keV(E, M_XE) - T_UMBRAL, 1.0, 20.0)
ax.axvline(E_cruce, color=C_XE, lw=1.0, ls=":", alpha=0.7)
ax.plot(E_cruce, T_UMBRAL, "o", ms=6, color=C_XE, zorder=5)
ax.annotate(rf"$E_\nu={coma(E_cruce,2)}$ MeV", (E_cruce, T_UMBRAL),
            xytext=(8.2, 2.25), fontsize=8.0, color=C_XE, ha="left", va="bottom",
            arrowprops=dict(arrowstyle="-", lw=0.6, color=C_XE))

ax.set_xlim(0, 12)
ax.set_ylim(0, None)
ax.set_xlabel(r"Energía del antineutrino, $E_\nu$  [MeV]")
ax.set_ylabel(r"$T^{\mathrm{máx}}(E_\nu)$  [keV]")
ax.legend(loc="upper left", fontsize=8.0)
formatear_coma(ax)
fig.tight_layout()
out = f"{BASE}/fig_tmax_kinematica.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
dest = f"{DEST}/fig_tmax_kinematica.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"  -> {dest}  (md5 {md5(out)})")

print("=" * 60)
print(f" M_Xe = {M_XE:.2f} MeV   M_Ar = {M_AR:.2f} MeV")
for E in (3.0, 5.0, 8.0):
    print(f"  E_nu={E:.0f} MeV: T_max Xe={Tmax_keV(E, M_XE):.4f} keV   "
          f"T_max Ar={Tmax_keV(E, M_AR):.4f} keV   "
          f"razon Ar/Xe={Tmax_keV(E, M_AR)/Tmax_keV(E, M_XE):.3f}")
print(f"\n  E_nu donde T_max^Xe cruza {T_UMBRAL} keV (umbral ~4 e-): {E_cruce:.4f} MeV")
print(f"  {out}")
