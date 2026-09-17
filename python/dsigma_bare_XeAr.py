#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dsigma/dT "desnuda" (F^2=1, sin flujo) del SM, Xe y Ar, a 3 energias de
neutrino fijas (Eν = 3, 5, 8 MeV), cada curva solo hasta su propio T_max(Eν).

Aisla la escala global de la seccion eficaz (~Q_W^2 M) y su forma casi plana
en T, antes de que entre el flujo o cualquier efecto de detector.

Motor real: se compila un driver minimo contra los modulos SIN modificarlos
(constants.f90 + xsections_nest.f90, funcion dsigma_dT, USE_HELM sin definir
-> F^2=1 por defecto), uno por carpeta (Xe: N_EventosCEvNS_NSIXe, Ar:
N_EventosCEvNS_NSIAr). Mismo Q_W que usa chi2_ideal_nest.f90
(QW_SM = -N/2 + (1-4 sin^2th_W)/2 * Z; N de Xe = <A>-Z, mezcla natural).

Salida: datos/dsigma_bare_{Xe,Ar}.dat, datos/fig_dsigma_dT_bare.png
"""
import os
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR
aplicar()

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DIR = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe", "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}
ENU_LIST = [3.0, 5.0, 8.0]   # MeV
N_T = 300

DRIVER = r"""
program dump_bare
  use constants, only: dp, N_Ge, Z_Ge, s2w, QV2, M_Ge
  use xsections_nest, only: dsigma_dT
  implicit none
  integer :: u, ie, it, n_t
  real(dp) :: QW_SM, E_nu, Tmax, T, dsdT
  real(dp), parameter :: enus(3) = (/ 3.0_dp, 5.0_dp, 8.0_dp /)

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2   = QW_SM**2
  n_t = __NT__

  open(newunit=u, file='__OUT__', status='replace')
  write(u,'(A)') '# E_nu[MeV]  T[keV]  dsigma_dT[cm^2/keV]  (F^2=1)'
  do ie = 1, 3
     E_nu = enus(ie)
     Tmax = 2.0_dp*E_nu**2 / (M_Ge + 2.0_dp*E_nu)     ! MeV, formula exacta
     do it = 1, n_t
        T = Tmax * real(it, dp) / real(n_t, dp)         ! evita T=0 y T=Tmax exactos
        dsdT = dsigma_dT(E_nu, T)
        write(u,'(F6.2, ES14.6, ES16.6)') E_nu, T*1000.0_dp, dsdT
     end do
  end do
  close(u)
end program dump_bare
"""


def run_target(tag):
    d = DIR[tag]
    out_dat = f"{BASE}/dsigma_bare_{tag}.dat"
    with tempfile.TemporaryDirectory() as tmp:
        src = os.path.join(tmp, "dump_bare.f90")
        exe = os.path.join(tmp, "dump_bare")
        code = DRIVER.replace("__OUT__", out_dat).replace("__NT__", str(N_T))
        with open(src, "w") as f:
            f.write(code)
        subprocess.run(
            ["gfortran", "-O2", "-ffree-line-length-none", "-o", exe,
             f"{d}/constants.f90", f"{d}/xsections_nest.f90", src],
            check=True, cwd=tmp)
        subprocess.run([exe], check=True, cwd=tmp)
    print(f"  -> {out_dat}")
    return out_dat


files = {tag: run_target(tag) for tag in ("Xe", "Ar")}

fig, ax = plt.subplots(figsize=(7.2, 5.2))
LS = {3.0: ":", 5.0: "--", 8.0: "-"}
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    d = np.loadtxt(files[tag], comments="#")
    for enu in ENU_LIST:
        sel = (d[:, 0] == enu) & (d[:, 2] > 3e-41)   # recorta el borde T->Tmax
        T, dsdT = d[sel, 1], d[sel, 2]                # (dsigma/dT->0, no aporta en log)
        order = np.argsort(T)
        ax.plot(T[order], dsdT[order], color=col, ls=LS[enu], lw=1.7)

# leyenda compuesta: color = blanco, estilo de linea = energia
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    ax.plot([], [], color=col, lw=1.7, label=tag)
for enu in ENU_LIST:
    ax.plot([], [], color="0.35", ls=LS[enu], lw=1.4, label=fr"$E_\nu={enu:.0f}$ MeV")

ax.set_yscale("log")
ax.set_ylim(3e-41, 4e-39)
ax.set_xlabel(r"Energía de retroceso nuclear, $T$  [keV]")
ax.set_ylabel(r"$d\sigma/dT$  [cm$^2$/keV]  ($F^2=1$)")
ax.legend(loc="upper right", ncol=2, fontsize=8.0)
fig.tight_layout()
out = f"{BASE}/fig_dsigma_dT_bare.png"
fig.savefig(out)
print(f"\n  {out}")
