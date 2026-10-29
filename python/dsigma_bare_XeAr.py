#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dsigma/dT "desnuda" (sin flujo) del SM, Xe y Ar, a 3 energias de neutrino
fijas (E_nu = 3, 5, 8 MeV), cada curva solo hasta su propio T_max(E_nu).
F^2=1 (continuo) y factor de forma de Helm (discontinuo, misma funcion
helm_F2 de xsections_nest.f90, llamada directamente -- sin depender de la
variable de entorno USE_HELM). Aisla la escala global de la seccion eficaz
(~Q_W^2 M) y su forma casi plana en T antes de que entre el flujo o
cualquier efecto de detector.

Panel principal: dsigma/dT, log-y, con una linea vertical punteada en cada
T_max. Subpanel: cociente Helm/F^2=1 (=helm_F2(T)) para las mismas curvas.

Motor real: se compila un driver minimo contra constants.f90 +
xsections_nest.f90 SIN modificarlos, uno por carpeta. Mismo Q_W que usa
chi2_ideal_nest.f90.

Salida: datos/dsigma_bare_{Xe,Ar}.dat, datos/fig_dsigma_dT_bare.png (+ copia en figuras_metodologia)
"""
import os
import shutil
import hashlib
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, formatear_coma
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0, "lines.linewidth": 1.3,
})

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
DIR = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe", "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}
ENU_LIST = [3.0, 5.0, 8.0]   # MeV
N_T = 300

DRIVER = r"""
program dump_bare
  use constants, only: dp, N_Ge, Z_Ge, s2w, QV2, M_Ge
  use xsections_nest, only: dsigma_dT, helm_F2
  implicit none
  integer :: u, ie, it, n_t
  real(dp) :: QW_SM, E_nu, Tmax, T, dsdT, f2
  real(dp), parameter :: enus(3) = (/ 3.0_dp, 5.0_dp, 8.0_dp /)

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2   = QW_SM**2
  n_t = __NT__

  open(newunit=u, file='__OUT__', status='replace')
  write(u,'(A)') '# E_nu[MeV]  T[keV]  dsigma_dT_F1[cm^2/keV]  dsigma_dT_Helm[cm^2/keV]  Tmax[keV]'
  do ie = 1, 3
     E_nu = enus(ie)
     Tmax = 2.0_dp*E_nu**2 / (M_Ge + 2.0_dp*E_nu)     ! MeV, formula exacta
     do it = 1, n_t
        T = Tmax * real(it, dp) / real(n_t, dp)         ! evita T=0 y T=Tmax exactos
        dsdT = dsigma_dT(E_nu, T)                        ! F^2=1 (USE_HELM no definida)
        f2 = helm_F2(T)
        write(u,'(F6.2, ES14.6, 2ES16.6, ES14.6)') E_nu, T*1000.0_dp, dsdT, dsdT*f2, Tmax*1000.0_dp
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

fig, (ax, axr) = plt.subplots(2, 1, figsize=(4.60, 4.3), sharex=False,
                                gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.32))
LS = {3.0: ":", 5.0: "--", 8.0: "-"}
tmax8 = {}
F2_tmax8 = {}
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    d = np.loadtxt(files[tag], comments="#")
    for enu in ENU_LIST:
        sel = (d[:, 0] == enu)   # sin piso: se ve caer a 0 en T_max, ylim recorta visualmente
        T, dsF1, dsHelm = d[sel, 1], d[sel, 2], d[sel, 3]
        order = np.argsort(T)
        T, dsF1, dsHelm = T[order], dsF1[order], dsHelm[order]
        ax.plot(T, dsF1, color=col, ls=LS[enu], lw=1.8)   # solo F^2=1 en el panel principal
        axr.plot(T, dsHelm / dsF1, color=col, ls=LS[enu], lw=1.5)   # F^2(T) de Helm; no depende de E_nu
        Tmax_keV = d[sel, 4][0]
        ax.axvline(Tmax_keV, color=col, lw=0.6, ls=":", alpha=0.4)
        if enu == 8.0:
            tmax8[tag] = Tmax_keV
            F2_tmax8[tag] = (dsHelm / dsF1)[-2]   # F2(T) es suave en T_max; evita dividir por dsF1->0 en el ultimo punto

for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    ax.plot([], [], color=col, lw=1.8, label=tag)
for enu in ENU_LIST:
    ax.plot([], [], color="0.35", ls=LS[enu], lw=1.4, label=rf"$E_\nu={enu:.0f}$ MeV")

ax.set_yscale("log")
ax.set_ylim(1e-42, 4e-39)
ax.set_ylabel(r"$d\sigma/dT$  [cm$^2$/keV]  ($F^2=1$)")
ax.legend(loc="upper right", ncol=2, fontsize=8.0)
formatear_coma(ax, y=False)

axr.set_ylabel(r"$F^2(T)$ de Helm")
axr.set_xlabel(r"Energía de retroceso nuclear, $T$  [keV]")
axr.axhline(1.0, color="0.7", lw=0.8)
axr.plot([], [], color="none", label="(no depende de $E_\\nu$)")
axr.legend(loc="upper right", fontsize=8.0, handlelength=0, handletextpad=0)
formatear_coma(axr)

fig.tight_layout()
out = f"{BASE}/fig_dsigma_dT_bare.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
dest = f"{DEST}/fig_dsigma_dT_bare.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"\n  T_max(8 MeV): Xe={tmax8['Xe']:.2f} keV  Ar={tmax8['Ar']:.2f} keV")
print(f"  F^2 de Helm en T_max(8 MeV):  Xe={F2_tmax8['Xe']:.3f}   Ar={F2_tmax8['Ar']:.3f}")
print(f"  -> {dest}  (md5 {md5(out)})")
