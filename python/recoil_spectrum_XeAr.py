#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Espectro de retroceso nuclear CEvNS del modelo estandar, dR/dT_nr, para Xe y
Ar (flujo hibrido Kopeikin U Huber-Mueller), CON el factor de forma de Helm
(el default de produccion; el efecto de F^2 queda en las figuras de seccion
eficaz, no aqui). Es el punto de partida del pipeline: lo que se integra
contra la seccion eficaz antes de NEST.

Panel inferior: cociente (dR/dT)_Xe / (dR/dT)_Ar, que tiende a
Q_W(Xe)^2/Q_W(Ar)^2 = 12,5 cuando T->0 (limite coherente, point-like) y
cruza 1 en el T donde la coherencia de Xe deja de compensar que sus
retrocesos son mas blandos (cinematica, Ec. T_max(E_nu)).

Motor: driver Fortran minimo compilado al vuelo contra constants.f90 +
xsections_nest.f90 + flux.f90 SIN modificarlos (mismo patron que
dsigma_bare_XeAr.py), integrando sobre el flujo real para cada T_nr.

Salida: datos/dR_dT_helm_{Xe,Ar}.dat, datos/fig_recoil_XeAr.png
        (+ copia en figuras_metodologia)
"""
import os
import shutil
import hashlib
import subprocess
import tempfile
from scipy.optimize import brentq
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, GRIS, NEGRO, formatear_coma, coma
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8, "lines.linewidth": 1.15,
})
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DIR = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe", "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}
T_PISO_NEST = 0.2    # keV, piso del modelo de ionizacion (NEST/LArNEST)
T_MIN_KEV = 0.01     # keV, piso de integracion/graficado (dR/dT no depende de NEST)

DRIVER = r"""
program recoil_helm
  use constants
  use xsections_nest, only: dsigma_dT, helm_F2
  use flux, only: flujo_diferencial, E_nu_max
  implicit none
  integer, parameter :: n_T = 800, n_E = 3000
  real(dp) :: T_nr, T_nr_min, T_nr_max, logT_min, logT_max, dlogT, E_nu, dE, peso, dummy
  real(dp) :: rate_F1, rate_Helm, atoms_per_kg, sec_per_day, f2, QW_SM
  integer :: i_T, i_E, u

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2
  dummy = flujo_diferencial(1.0_dp)
  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp

  T_nr_min = __TMIN__ / 1000.0_dp
  T_nr_max = 2.0_dp*E_nu_max**2/(M_Ge+2.0_dp*E_nu_max)
  logT_min = log(T_nr_min); logT_max = log(T_nr_max)
  dlogT = (logT_max - logT_min) / real(n_T-1, dp)
  dE = E_nu_max / real(n_E-1, dp)

  open(newunit=u, file='__OUT__', status='replace')
  write(u,'(A)') '# T_nr[keV]  dRdT_F1[ev/(kg dia keV)]  dRdT_Helm[ev/(kg dia keV)]'
  do i_T = 1, n_T
     T_nr = exp(logT_min + (i_T-1)*dlogT)
     f2 = helm_F2(T_nr)
     rate_F1 = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E-1)*dE
        if (E_nu < sqrt(M_Ge*T_nr/2.0_dp)) cycle
        peso = merge(0.5_dp, 1.0_dp, i_E==1 .or. i_E==n_E)
        rate_F1 = rate_F1 + flujo_diferencial(E_nu)*dsigma_dT(E_nu,T_nr)*dE*peso
     end do
     rate_F1  = rate_F1 * atoms_per_kg * sec_per_day
     rate_Helm = rate_F1 * f2
     write(u,'(F9.5, 2ES16.6)') T_nr*1000.0_dp, rate_F1, rate_Helm
  end do
  close(u)
end program recoil_helm
"""


def run_target(tag):
    d = DIR[tag]
    out_dat = f"{BASE}/dR_dT_helm_{tag}.dat"
    with tempfile.TemporaryDirectory() as tmp:
        src, exe = os.path.join(tmp, "recoil_helm.f90"), os.path.join(tmp, "recoil_helm")
        code = DRIVER.replace("__OUT__", out_dat).replace("__TMIN__", str(T_MIN_KEV))
        with open(src, "w") as f:
            f.write(code)
        subprocess.run(["gfortran", "-O2", "-ffree-line-length-none", "-o", exe,
                         f"{d}/constants.f90", f"{d}/xsections_nest.f90", f"{d}/flux.f90", src],
                        check=True, cwd=tmp)
        subprocess.run([exe], check=True, cwd=tmp)
    print(f"  -> {out_dat}")
    return out_dat


files = {tag: run_target(tag) for tag in ("Xe", "Ar")}
dat = {tag: np.loadtxt(files[tag], comments="#") for tag in ("Xe", "Ar")}

T_MIN_PLOT, T_MAX_PLOT = 0.05, 5.0

fig, (ax, axr) = plt.subplots(2, 1, figsize=(4.61, 4.4),
                                gridspec_kw=dict(height_ratios=[2.0, 1.0], hspace=0.10))
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    T, rHelm = dat[tag][:, 0], dat[tag][:, 2]
    m = rHelm > 1e-6
    ax.plot(T[m], rHelm[m], color=col, lw=1.6, label=tag)

# cociente Xe/Ar (interpolado a una malla log comun)
Tg = np.logspace(np.log10(T_MIN_PLOT), np.log10(T_MAX_PLOT), 600)
rXe = np.interp(Tg, dat["Xe"][:, 0], dat["Xe"][:, 2], left=np.nan, right=np.nan)
rAr = np.interp(Tg, dat["Ar"][:, 0], dat["Ar"][:, 2], left=np.nan, right=np.nan)
ratio = rXe / rAr
QW2_ref = 12.516   # (Q_W,Xe/Q_W,Ar)^2, point-like, T->0
m2 = np.isfinite(ratio) & (ratio > 0)
axr.plot(Tg[m2], ratio[m2], color=NEGRO, lw=1.4, label=r"$(dR/dT)_{\rm Xe}/(dR/dT)_{\rm Ar}$")
axr.axhline(QW2_ref, color=GRIS, lw=0.8, ls=":")
axr.text(T_PISO_NEST * 1.25, QW2_ref * 1.08, rf"$Q_W^2({{\rm Xe}})/Q_W^2({{\rm Ar}})={coma(QW2_ref,1)}$",
         fontsize=6.8, color="0.35", va="bottom")
axr.axhline(1.0, color=GRIS, lw=0.8, ls=":")

f = lambda logT: np.interp(10**logT, Tg[m2], ratio[m2]) - 1.0
T_cross = 10 ** brentq(f, np.log10(T_MIN_PLOT), np.log10(1.0))
axr.plot(T_cross, 1.0, "o", ms=5, mfc="white", mec=NEGRO, mew=1.3, zorder=6)
axr.annotate(rf"$T={coma(T_cross,2)}$ keV", (T_cross, 1.0), textcoords="offset points",
             xytext=(8, 6), fontsize=6.8, color=NEGRO)

for a in (ax, axr):
    a.axvspan(T_MIN_PLOT, T_PISO_NEST, color="0.85", alpha=0.6, zorder=0)
    a.axvline(T_PISO_NEST, color="0.4", lw=0.9, ls=":")
    a.set_xscale("log")
    a.set_xlim(T_MIN_PLOT, T_MAX_PLOT)

ax.text(np.sqrt(T_MIN_PLOT * T_PISO_NEST), 0.55, "por debajo del\nmínimo de NEST",
        transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7.0, color="0.35")
ax.set_yscale("log")
ax.set_ylabel(r"$dR/dT$  [eventos / (kg$\cdot$día$\cdot$keV)]")
ax.legend(loc="upper right", fontsize=7.8)
plt.setp(ax.get_xticklabels(), visible=False)
ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterSciNotation())
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

axr.set_yscale("log")
axr.set_ylim(0.3, 28)
axr.set_ylabel(r"$(dR/dT)_{\rm Xe}\,/\,(dR/dT)_{\rm Ar}$", fontsize=8.2)
axr.set_xlabel(r"Energía de retroceso nuclear  $T$  [keV]")

formatear_coma(axr)

fig.tight_layout()
out = f"{BASE}/fig_recoil_XeAr.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = f"{DEST}/fig_recoil_XeAr.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"  T donde (dR/dT)_Xe/(dR/dT)_Ar cruza 1: {T_cross:.4f} keV")
i02 = np.argmin(np.abs(Tg - 0.2))
print(f"  cociente en T=0,2 keV: {ratio[i02]:.4f}")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
