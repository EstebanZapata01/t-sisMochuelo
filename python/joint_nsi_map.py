#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapa NSI conjunto (idea 4): CONUS+ (real) (+) Xe (real), con Ar (proyeccion
Asimov) mostrado APARTE y etiquetado, siguiendo la recomendacion ya escrita
en doc/metodologia.tex Sec. "Combinacion de experimentos: real (+) real vs.
real (+) Asimov": el mapa honesto de "estado actual" es CONUS+ (+) Xe-real
(dos experimentos con datos observados, suma de chi2 = -2 ln del producto de
verosimilitudes); Ar solo se agrega como capa etiquetada de proyeccion, no
sin asterisco.

Construido sobre las grillas que produce chi2_nsi_generic (ver
python/../FORTRAN90/chi2_nsi_generic/), ya verificadas bit-a-bit contra
chi2_nsi_2DXe.dat y chi2_nsi_2Dconus.dat (diff maximo 1e-5, precision de
formato del .dat) -- sumarlas punto a punto es valido porque comparten
exactamente la misma malla 1000x1000 en [-1,1]^2 y el mismo ipar=5.

Entradas (datos/):
    generic_Xe_nsi2D.dat       chi2_Xe,real(eps)     (de chi2_nsi_generic)
    generic_conus_nsi2D.dat    chi2_CONUS+(eps)      (de chi2_nsi_generic)
    chi2_nsi_2DAr_ideal.dat    Delta_chi2_Ar,Asimov(eps)  (de chi2_ideal_nest.f90)

Salidas (datos/):
    chi2_nsi_conjunto.dat      eps_x eps_y chi2_actual chi2_actual+Ar_proyectado
    fig_nsi_conjunto.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9, "axes.grid": True,
    "grid.color": "0.82", "grid.linewidth": 0.5,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "legend.frameon": False, "legend.fontsize": 9,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})

def load_labels(path):
    with open(path, encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    return lines[0], lines[1]

xlabel, ylabel = load_labels(f"{BASE}/nsi_config_conus.txt")

# --------------------------------------------------------------- cargar grillas
xe = np.loadtxt(f"{BASE}/generic_Xe_nsi2D.dat", comments="#")
co = np.loadtxt(f"{BASE}/generic_conus_nsi2D.dat", comments="#")
ar = np.loadtxt(f"{BASE}/chi2_nsi_2DAr_ideal.dat", comments="#")

assert np.allclose(xe[:, 0], co[:, 0]) and np.allclose(xe[:, 1], co[:, 1]), \
    "las grillas Xe/CONUS+ no comparten malla -- no se pueden sumar punto a punto"
assert np.allclose(xe[:, 0], ar[:, 0]) and np.allclose(xe[:, 1], ar[:, 1]), \
    "la grilla de Ar no coincide con la de Xe/CONUS+"

eps_x_flat, eps_y_flat = xe[:, 0], xe[:, 1]
chi2_xe, chi2_co, dchi2_ar = xe[:, 2], co[:, 2], ar[:, 2]

n = int(round(np.sqrt(len(eps_x_flat))))
EX = eps_x_flat.reshape(n, n)
EY = eps_y_flat.reshape(n, n)
CHI2_XE = chi2_xe.reshape(n, n)
CHI2_CO = chi2_co.reshape(n, n)
DCHI2_AR = dchi2_ar.reshape(n, n)

# ---------------------------------------------------------- estado actual: real+real
chi2_actual = CHI2_CO + CHI2_XE
chi2_actual_min = chi2_actual.min()
dchi2_actual = chi2_actual - chi2_actual_min

# CONUS+ solo, para comparar cuanto aporta Xe realmente
chi2_co_min = CHI2_CO.min()
dchi2_co = CHI2_CO - chi2_co_min

CL90 = 4.605  # 2 g.d.l.
frac_actual = np.mean(dchi2_actual < CL90)
frac_co     = np.mean(dchi2_co < CL90)

print("=" * 78)
print(" MAPA NSI CONJUNTO: CONUS+(real) + Xe(real), Ar(Asimov) aparte")
print("=" * 78)
print(f" chi2_min CONUS+ solo      = {chi2_co_min:.4f}")
print(f" chi2_min CONUS+ + Xe,real = {chi2_actual_min:.4f}  "
      f"(chi2_min(Xe,real) solo = {CHI2_XE.min():.4f})")
print(f" Fraccion de |eps|<=1 con Dchi2<4.605 (90%,2gdl):")
print(f"   CONUS+ solo         = {frac_co*100:6.2f} %")
print(f"   CONUS+ + Xe,real    = {frac_actual*100:6.2f} %   "
      f"(Xe real aporta muy poco: A_90~107xSM no toca la caja |eps|<=1)")

# capa de Ar proyectado (Asimov): se suma SOLO para mostrar "cuanto apretaria"
# un futuro run de UAr, en un panel APARTE (no en el mapa "actual" sin asterisco)
chi2_actual_mas_ar = dchi2_actual + DCHI2_AR  # ambas ya "Delta" respecto a su propio minimo local
frac_mas_ar = np.mean(chi2_actual_mas_ar < CL90)
print(f"   CONUS+ + Xe,real + Ar,PROYECTADO (condicionada) = {frac_mas_ar*100:6.2f} %")

with open(f"{BASE}/chi2_nsi_conjunto.dat", "w", encoding="utf-8") as f:
    f.write("# eps_x  eps_y  Dchi2_actual(CONUS+ + Xe,real)  "
            "Dchi2_actual+Ar,proyectado(condicionada)\n")
    for i in range(n):
        for j in range(n):
            f.write(f"{EX[i,j]: .6E} {EY[i,j]: .6E} "
                     f"{dchi2_actual[i,j]: .6E} {chi2_actual_mas_ar[i,j]: .6E}\n")
        f.write("\n")

# --------------------------------------------------------------- figura
fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.4), sharey=True)

ax = axes[0]
ax.contourf(EX, EY, dchi2_actual, levels=[0, CL90, dchi2_actual.max()],
            colors=["0.55", "0.92"])
ax.contour(EX, EY, dchi2_actual, levels=[CL90], colors="black", linewidths=1.4)
ax.plot(0, 0, "+", color="black", ms=10, mew=1.6)
ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
ax.set_title(f"Estado actual: CONUS+ $\\oplus$ Xe (datos reales)\n"
             f"{frac_actual*100:.1f}% de $|\\varepsilon|\\leq 1$ permitido (90%)")
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1)

ax2 = axes[1]
ax2.contourf(EX, EY, dchi2_actual, levels=[0, CL90, dchi2_actual.max()],
             colors=["0.55", "0.92"])
ax2.contour(EX, EY, dchi2_actual, levels=[CL90], colors="black", linewidths=1.4)
ax2.contour(EX, EY, chi2_actual_mas_ar, levels=[CL90], colors="0.35",
            linewidths=1.6, linestyles="--")
ax2.plot(0, 0, "+", color="black", ms=10, mew=1.6)
ax2.set_xlabel(xlabel)
ax2.set_title("+ Ar proyectado (Asimov, condicionada)")
ax2.set_xlim(-1, 1); ax2.set_ylim(-1, 1)

legend_elems = [
    Patch(facecolor="0.55", edgecolor="black", label=r"Permitido (datos reales), $\Delta\chi^2<4.605$"),
    Line2D([0], [0], color="0.35", lw=1.6, ls="--",
           label="+ proyección Ar (Asimov) — condicionada, no medición"),
    Line2D([0], [0], marker="+", color="black", lw=0, ms=9, mew=1.6, label="SM ($\\varepsilon=0$)"),
]
fig.legend(handles=legend_elems, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Mapa NSI conjunto (plano $\\varepsilon_{ee}^{dV}$–$\\varepsilon_{e\\mu}^{dV}$, ipar=5)")
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(f"{BASE}/fig_nsi_conjunto.png")
print(f"\n  {BASE}/fig_nsi_conjunto.png")
print(f"  {BASE}/chi2_nsi_conjunto.dat")
