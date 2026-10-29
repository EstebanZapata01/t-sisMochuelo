#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cascada (waterfall) de por que el argon le gana al xenon.

La Tabla de cargas debiles dice que Xe tiene ventaja de coherencia
(x3.8 por unidad de masa); la comparacion ideal dice que Ar termina
ganando por hasta x254 con umbral N_e>=4. Este script muestra en que
paso exacto de la cadena se invierte la ventaja, panel por panel, todo
con datos que el pipeline ideal ya escribe (ningun calculo nuevo):

  (1) dR/dT           -- espectro de retroceso nuclear (teoria pura,
                          antes de NEST/LArNEST y de la extraccion).
                          Aqui manda la coherencia: Xe por encima.
  (2) R_k "creados"   -- tras el modelo binomial de Fano (yield +
                          fluctuacion de ionizacion). Xe cae en picado
                          hacia N_e alto porque su yield en el umbral
                          es ~0 y sus retrocesos son blandos.
  (3) R_k "extraidos" -- tras la extraccion de la interfase (EEE).
                          Xe (EEE ~ 0.33 efectiva sobre p_F) pierde casi
                          todo; Ar (EEE ~ 0.99) no pierde nada.
  (4) R_tot(N_e >= k) -- integrado. El punto final: la ventaja Ar/Xe.

Entradas (datos/):
  espectro_continuo{Xe,Ar}.dat   : T_nr[keV]  Kop  Mue  Comb
  ionization_electrones{,_Ar}.dat : N_e  creados  extraidos  [ev/(kg dia)]
  sensib_ideal_{Xe,Ar}.dat        : NE_LO F_mode R_tot ... (fila Fnest)

Salidas: datos/fig_waterfall_{1_retroceso,2_extraidos,3_integrado}.png
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
C_XE, C_AR = "#33546e", "#a86a43"



def load(path):
    return np.loadtxt(path, comments="#")


dRdT_xe = load(f"{BASE}/espectro_continuoXe.dat")   # T Kop Mue Comb
dRdT_ar = load(f"{BASE}/espectro_continuoAr.dat")
ne_xe = load(f"{BASE}/ionization_electrones.dat")    # N_e creados extraidos
ne_ar = load(f"{BASE}/ionization_electrones_Ar.dat")


def rtot_ge(path):
    """R_tot(N_e>=k) para k=1..4 de la fila Fnest de sensib_ideal_*."""
    k, r = [], []
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            if p[1] == "Fnest":
                k.append(int(p[0]))
                r.append(float(p[2]))
    return np.array(k), np.array(r)


k_xe, rtot_xe = rtot_ge(f"{BASE}/sensib_ideal_Xe.dat")
k_ar, rtot_ar = rtot_ge(f"{BASE}/sensib_ideal_Ar.dat")

from estilo_tesis import FIG15
aplicar(grande=True)

def nueva():
    return plt.subplots(figsize=FIG15)

def guarda(fig, nombre, copiar=False):
    fig.tight_layout()
    out = f"{BASE}/{nombre}.png"
    fig.savefig(out, dpi=300); plt.close(fig)
    print(f"  -> {out}")
    if copiar:
        dest = f"{DEST}/{nombre}.png"
        shutil.copyfile(out, dest)
        md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
        assert md5(out) == md5(dest)
        print(f"  -> {dest}  (md5 {md5(out)})")

# ---- (1): dR/dT, teoria pura, con la ROI sombreada
fig, a = nueva()
a.plot(dRdT_xe[:, 0], dRdT_xe[:, 3], color=C_XE, label="Xe")
a.plot(dRdT_ar[:, 0], dRdT_ar[:, 3], color=C_AR, label="Ar")
a.axvspan(0.2, 1.0, color="0.9", zorder=0)
a.set_yscale("log"); a.set_xlim(0, 4)
a.set_xlabel(r"$T_{\rm nr}$ [keV]")
a.set_ylabel(r"$dR/dT$ [ev/(kg$\cdot$día$\cdot$keV)]")
a.legend(loc="upper right")
a.text(0.6, 0.06, "ROI CE$\\nu$NS\n(0,2–1 keV)", transform=a.get_xaxis_transform(),
       ha="center", va="bottom", fontsize=10.5, color="0.35")
guarda(fig, "fig_waterfall_1_retroceso")

# ---- (2): extraidos (creados punteados)
fig, a = nueva()
a.step(ne_xe[:, 0], ne_xe[:, 2], where="mid", color=C_XE, label="Xe, extraídos")
a.step(ne_ar[:, 0], ne_ar[:, 2], where="mid", color=C_AR, label="Ar, extraídos")
a.step(ne_xe[:, 0], ne_xe[:, 1], where="mid", color=C_XE, lw=1.2, ls=":", alpha=0.7, label="Xe, creados")
a.step(ne_ar[:, 0], ne_ar[:, 1], where="mid", color=C_AR, lw=1.2, ls=":", alpha=0.7, label="Ar, creados")
a.set_yscale("log"); a.set_xlim(0, 15); a.set_ylim(1e-4, 5e1)
a.set_xlabel(r"$N_e$")
a.set_ylabel(r"$R_k$ [ev/(kg$\cdot$día)]")
a.legend(loc="upper right")
guarda(fig, "fig_waterfall_2_extraidos")

# ---- (3): 6 categorias -- "todos los retrocesos" (T:0->Tmax(E), con Helm),
#           "T>0,2 keV" (integral del archivo, mismo piso que la tabla NEST),
#           N_e>=1..4 -- con Xe(EEE=0.328) y Ar; subpanel: cociente Ar/Xe.
import os
import subprocess
import tempfile
from estilo_tesis import formatear_coma, coma, GRIS

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
DIRF = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe", "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}

DRIVER_TODOS = r"""
program todos_retrocesos
  use constants
  use flux, only: flujo_diferencial, E_nu_max
  use xsections_nest, only: dsigma_dT
  implicit none
  integer, parameter :: n_E = 2000, n_T = 2000
  real(dp) :: E_nu, dE, pesoE, QW_SM, Tmax, dT, T_MeV, sigma_tot, R_tot, dummy
  real(dp) :: atoms_per_kg, sec_per_day
  integer :: i_E, i_T
  dummy = flujo_diferencial(1.0_dp)
  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2
  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp
  dE = E_nu_max / real(n_E-1, dp)
  ! sigma_total(E) = integral real de dsigma_dT (Helm por defecto, USE_HELM=0
  ! la apaga), T:0->Tmax(E) exacto, NO la aproximacion analitica F=1.
  R_tot = 0.0_dp
  do i_E = 1, n_E
     E_nu = (i_E-1)*dE
     if (E_nu <= 0.0_dp) cycle
     Tmax = 2.0_dp * E_nu**2 / (M_Ge + 2.0_dp * E_nu)
     dT = Tmax / n_T
     sigma_tot = 0.0_dp
     do i_T = 1, n_T
        T_MeV = (real(i_T,dp) - 0.5_dp) * dT
        ! dsigma_dT da cm^2/keV; integrar en dT_keV, no en dT[MeV].
        sigma_tot = sigma_tot + dsigma_dT(E_nu, T_MeV) * dT * 1000.0_dp
     end do
     pesoE = merge(0.5_dp, 1.0_dp, i_E==1 .or. i_E==n_E)
     R_tot = R_tot + flujo_diferencial(E_nu) * sigma_tot * dE * pesoE
  end do
  R_tot = R_tot * atoms_per_kg * sec_per_day
  write(*,'(ES20.10)') R_tot
end program todos_retrocesos
"""


def todos_los_retrocesos(tag):
    d = DIRF[tag]
    with tempfile.TemporaryDirectory() as tmp:
        src, exe = os.path.join(tmp, "tr.f90"), os.path.join(tmp, "tr")
        with open(src, "w") as f:
            f.write(DRIVER_TODOS)
        subprocess.run(["gfortran", "-O2", "-ffree-line-length-none", "-o", exe,
                         f"{d}/constants.f90", f"{d}/flux.f90", f"{d}/xsections_nest.f90", src],
                        check=True, cwd=tmp)
        out = subprocess.run([exe], check=True, cwd=tmp, capture_output=True, text=True)
    return float(out.stdout.strip().splitlines()[-1])


R_todos_xe = todos_los_retrocesos("Xe")
R_todos_ar = todos_los_retrocesos("Ar")

# "T>0,2 keV" = integral completa de dR/dT desde el piso de la tabla NEST (espectro_continuo{Xe,Ar}.dat).
R_02_xe = np.trapz(dRdT_xe[:, 3], dRdT_xe[:, 0])
R_02_ar = np.trapz(dRdT_ar[:, 3], dRdT_ar[:, 0])

CAT = ["todos los\nretrocesos", "$T\\!>\\!0{,}2$\nkeV", r"$N_e\!\geq\!1$",
       r"$N_e\!\geq\!2$", r"$N_e\!\geq\!3$", r"$N_e\!\geq\!4$"]
x = np.arange(len(CAT))
Rxe_full = np.concatenate(([R_todos_xe, R_02_xe], rtot_xe))
Rar_full = np.concatenate(([R_todos_ar, R_02_ar], rtot_ar))

aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0, "lines.linewidth": 1.3,
})
fig, (a, ar_) = plt.subplots(2, 1, figsize=(5.22, 4.8),
                              gridspec_kw=dict(height_ratios=[2.3, 1.0], hspace=0.08))

# Banda del Qy del argon: LArNEST-ReD (linea, igual que siempre) hasta la
# envolvente de DarkSide-50 central +-1sigma (sin correlacion), solo en
# N_e>=1..4 (arXiv:2107.08087 Sec. V; ver ar_qy_darkside50.py).
try:
    rds = np.loadtxt(f"{BASE}/rtot_darkside50_Ar.dat", comments="#", usecols=(1, 2, 3, 4))
    ds_c, ds_p1, ds_m1 = rds[0], rds[1], rds[2]
    ds_lo = np.minimum(ds_m1, ds_p1)
    ds_hi = np.maximum(ds_m1, ds_p1)
    x_banda = x[2:]   # N_e>=1..4
    a.fill_between(x_banda, np.minimum(rtot_ar, ds_lo), np.maximum(rtot_ar, ds_hi),
                    color=C_AR, alpha=0.35, lw=0,
                    label=r"Ar, banda de $Q_y$ (LArNEST–ReD a DS-50 $\pm1\sigma$)")
    a.plot(x_banda, ds_c, "--", color=C_AR, lw=1.3, alpha=0.9, label="Ar, DS-50 (central)")
    tiene_banda_ds50 = True
except OSError:
    tiene_banda_ds50 = False

a.plot(x, Rxe_full, "o-", color=C_XE, label=r"Xe ($\mathrm{EEE}=0{,}328$)")
a.plot(x, Rar_full, "o-", color=C_AR, label=r"Ar ($\mathrm{EEE}=0{,}99$)")
a.set_yscale("log")
a.set_ylim(Rxe_full.min() * 0.3, Rar_full.max() * 6)
a.set_ylabel(r"$R_{\rm tot}$ [ev/(kg$\cdot$día)]")
a.legend(loc="lower left", fontsize=8.0)
plt.setp(a.get_xticklabels(), visible=False)

ratio_full = Rar_full / Rxe_full
if tiene_banda_ds50:
    ratio_lo = np.minimum(rtot_ar, ds_lo) / rtot_xe
    ratio_hi = np.maximum(rtot_ar, ds_hi) / rtot_xe
    ar_.fill_between(x_banda, ratio_lo, ratio_hi, color=C_AR, alpha=0.35, lw=0)
    ar_.plot(x_banda, ds_c / rtot_xe, "--", color=C_AR, lw=1.3, alpha=0.9)
ar_.plot(x, ratio_full, "o-", color="0.25", label="Ar / Xe (EEE=0,328)")
ar_.axhline(1.0, color=GRIS, lw=1.0, ls=":")
ar_.set_yscale("log")
ar_.set_ylim(ratio_full.min() * 0.22, ratio_full.max() * 12.0)
ar_.set_xticks(x); ar_.set_xticklabels(CAT, fontsize=8.0)
ar_.set_ylabel(r"Ar / Xe")
ar_.legend(loc="upper left", fontsize=8.0)
for xx, rr in zip(x, ratio_full):
    txt = rf"$\times{coma(rr,2)}$" if rr < 1 else (rf"$\times{coma(rr,1)}$" if rr < 10 else rf"$\times{rr:.0f}$")
    if xx == 0:
        ar_.annotate(txt, (xx, rr), textcoords="offset points", xytext=(14, -3),
                     ha="left", va="top", fontsize=8.0, color="0.2")
    else:
        ar_.annotate(txt, (xx, rr), textcoords="offset points", xytext=(0, 7),
                     ha="center", fontsize=8.0, color="0.2")

formatear_coma(ar_, x=False)
fig.tight_layout()
guarda(fig, "fig_waterfall_3_integrado", copiar=True)

print("=" * 74)
print(" CASCADA Xe vs Ar")
print("=" * 74)
print(f"  todos los retrocesos (T:0->Tmax(E), con Helm): Xe {R_todos_xe:.3f}  Ar {R_todos_ar:.3f}"
      f"  Xe/Ar x{R_todos_xe/R_todos_ar:.3f}  Ar/Xe x{R_todos_ar/R_todos_xe:.4f}")
print(f"  fraccion de la tasa total que sobrevive T>0,2 keV: Xe {R_02_xe/R_todos_xe:.4f}"
      f"   Ar {R_02_ar/R_todos_ar:.4f}")
print(f"  T>0,2 keV: Xe {R_02_xe:.3f}  Ar {R_02_ar:.3f}  Ar/Xe x{R_02_ar/R_02_xe:.3f}")
ratio = rtot_ar / rtot_xe
print(f"  R_tot(N_e>=1): Xe(0,328) {rtot_xe[0]:.2f}  Ar {rtot_ar[0]:.2f}   Ar/Xe x{ratio[0]:.2f}")
print(f"  R_tot(N_e>=4): Xe(0,328) {rtot_xe[3]:.3f}  Ar {rtot_ar[3]:.2f}   Ar/Xe x{ratio[3]:.0f}")
print(f"  Xe pierde en la extracción (N_e=4): creados {ne_xe[4,1]:.3f} -> extraídos {ne_xe[4,2]:.4f}")
print(f"  Ar conserva (N_e=4):                creados {ne_ar[4,1]:.3f} -> extraídos {ne_ar[4,2]:.3f}")
