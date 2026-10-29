#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Distribuciones normalizadas en E_nu: el flujo del reactor, los eventos de
Xe con N_e verdadero >= 1 (cadena ideal), los eventos de Xe con N_e
verdadero >= 4 (aproximacion de la ROI real) y los eventos de Ar con N_e
verdadero >= 1 (ROI 1-5, casi todo el rango util). Es decir,
Phi(E) * integral dT dsigma/dT * P(N_e>=k | T), normalizado a area 1.

Mensaje: los eventos del xenon en la ROI de RED-100 vienen de la cola alta
del espectro del reactor; el argon muestrea el grueso del espectro.

Motor real: driver Fortran compilado al vuelo contra constants.f90 +
xsections_nest.f90 + Tnr_to_e.f90 + mod_stats.f90 + flux.f90 SIN
modificarlos, uno por carpeta (mismo patron que mainred100_nest.f90, pero
con el orden de integracion invertido: E_nu afuera, T adentro, para
obtener dR/dE_nu en vez de dR/dT).

Salida: datos/energia_efectiva_{Xe,Ar}.dat, datos/fig_energia_efectiva.png
"""
import os
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, NEGRO, GRIS, coma, formatear_coma
aplicar(grande=True)

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DIR = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe", "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}

DRIVER = r"""
program energia_efectiva
  use constants
  use mod_tnr_to_e
  use xsections_nest, only: dsigma_dT
  use mod_stats, only: binomial_prob
  use flux, only: flujo_diferencial, E_nu_max
  implicit none
  integer, parameter :: n_E = 400, n_T = 200, n_ion = 100
  real(dp) :: E_nu, dE, T_min_keV, T_min, T_max, dT, dsdT, w, dummy
  real(dp) :: rate1, rate4, Psel1, Psel4, pk, flux_here
  integer :: i_E, i_T, k, n_F, u_out
  real(dp) :: p_F, atoms_per_kg, sec_per_day, QW_SM

  call inicializar_nest()
  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2
  dummy = flujo_diferencial(1.0_dp)

  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp
  T_min_keV = T_dense(1)
  T_min = T_min_keV / 1000.0_dp      ! MeV, mismas unidades que E_nu en dsigma_dT

  dE = E_nu_max / real(n_E - 1, dp)
  open(newunit=u_out, file='__OUT__', status='replace')
  write(u_out,'(A)') '# E_nu[MeV]  flujo[nu/cm2/s/MeV]  rate_Ne_ge1[ev/(kg dia MeV)]  rate_Ne_ge4[ev/(kg dia MeV)]'
  do i_E = 1, n_E
     E_nu = (i_E - 1) * dE
     flux_here = flujo_diferencial(E_nu)
     T_max = 2.0_dp * E_nu**2 / (M_Ge + 2.0_dp * E_nu)     ! MeV
     rate1 = 0.0_dp; rate4 = 0.0_dp
     if (T_max > T_min) then
        dT = (T_max - T_min) / real(n_T - 1, dp)
        do i_T = 1, n_T
           w = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)
           dsdT = dsigma_dT(E_nu, T_min + (i_T - 1) * dT)
           call obtener_nest_binomial((T_min + (i_T - 1) * dT) * 1000.0_dp, n_F, p_F)
           Psel1 = 0.0_dp; Psel4 = 0.0_dp
           do k = 1, min(n_ion, n_F)
              pk = binomial_prob(k, n_F, p_F * EEE)
              Psel1 = Psel1 + pk
              if (k >= 4) Psel4 = Psel4 + pk
           end do
           rate1 = rate1 + dsdT * Psel1 * dT * w
           rate4 = rate4 + dsdT * Psel4 * dT * w
        end do
     end if
     rate1 = rate1 * flux_here * atoms_per_kg * sec_per_day
     rate4 = rate4 * flux_here * atoms_per_kg * sec_per_day
     write(u_out,'(F9.5, 3ES16.6)') E_nu, flux_here, rate1, rate4
  end do
  close(u_out)
end program energia_efectiva
"""


def run_target(tag):
    d = DIR[tag]
    out_dat = f"{BASE}/energia_efectiva_{tag}.dat"
    with tempfile.TemporaryDirectory() as tmp:
        src = os.path.join(tmp, "energia_efectiva.f90")
        exe = os.path.join(tmp, "energia_efectiva")
        code = DRIVER.replace("__OUT__", out_dat)
        with open(src, "w") as f:
            f.write(code)
        subprocess.run(
            ["gfortran", "-O2", "-ffree-line-length-none", "-o", exe,
             f"{d}/constants.f90", f"{d}/mod_stats.f90", f"{d}/Tnr_to_e.f90",
             f"{d}/xsections_nest.f90", f"{d}/flux.f90", src],
            check=True, cwd=tmp)
        subprocess.run([exe], check=True, cwd=d)
    print(f"  -> {out_dat}")
    return out_dat


files = {tag: run_target(tag) for tag in ("Xe", "Ar")}

d_xe = np.loadtxt(files["Xe"], comments="#")
d_ar = np.loadtxt(files["Ar"], comments="#")
E_xe, flujo, xe_ne1, xe_ne4 = d_xe[:, 0], d_xe[:, 1], d_xe[:, 2], d_xe[:, 3]
E_ar, ar_ne1 = d_ar[:, 0], d_ar[:, 2]


def normalizar(E, y):
    area = np.trapz(y, E)
    return y / area


CURVAS = [
    ("Flujo del reactor", E_xe, normalizar(E_xe, flujo), NEGRO, "-"),
    (r"Xe, $N_e\geq1$ (cadena ideal)", E_xe, normalizar(E_xe, xe_ne1), C_XE, "-"),
    (r"Xe, $N_e\geq4$ (umbral de RED-100)", E_xe, normalizar(E_xe, xe_ne4), C_XE, "--"),
    (r"Ar, $N_e\geq1$", E_ar, normalizar(E_ar, ar_ne1), C_AR, "-"),
]


def mediana_y_cola(E, y):
    cdf = np.concatenate(([0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(E))))
    cdf /= cdf[-1]
    med = np.interp(0.5, cdf, E)
    frac8 = 1.0 - np.interp(8.0, E, cdf)
    return med, frac8


fig, ax = plt.subplots(figsize=(7.6, 6.0))
ax.axvspan(8.0, 10.0, color=GRIS, alpha=0.15, zorder=0)
ax.text(9.0, 0.55, "forma\nextrapolada", transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=8.6, color="0.35")

for lab, E, y, col, ls in CURVAS:
    ax.plot(E, y, color=col, ls=ls, lw=2.0, label=lab)

ax.set_xlim(0, 10)
ax.set_ylim(bottom=0)
ymax = ax.get_ylim()[1]
for lab, E, y, col, ls in CURVAS:
    med, frac8 = mediana_y_cola(E, y)
    y_med = np.interp(med, E, y)
    ax.plot([med, med], [y_med - 0.03 * ymax, y_med + 0.03 * ymax], color=col, lw=2.2, solid_capstyle="butt")

ax.set_xlabel(r"Energía del antineutrino, $E_\nu$  [MeV]")
ax.set_ylabel(r"Distribución normalizada  [MeV$^{-1}$]")
ax.legend(loc="upper right", fontsize=10)
formatear_coma(ax)

fig.tight_layout()
out = f"{BASE}/fig_energia_efectiva.png"
fig.savefig(out, dpi=300, bbox_inches="tight")

print("\n  Mediana de E_nu y fraccion con E_nu>8 MeV:")
for lab, E, y, col, ls in CURVAS:
    med, frac8 = mediana_y_cola(E, y)
    print(f"    {lab}: mediana E_nu = {med:.3f} MeV   fraccion(E_nu>8 MeV) = {frac8:.4f}")

print(f"  -> {out}")
