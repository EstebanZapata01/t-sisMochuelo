#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perfil Delta-chi2(A) construido por el MISMO binario (chi2_nsi_generic) para
Xe y para CONUS+, en un solo panel comparativo -- la prueba visual de la
validacion cruzada (Sec. "Validacion cruzada", doc/metodologia.tex): el
panel de Xe debe verse igual que fig9_chi2_perfil.pdf (construida por
chi2.f90 directamente); el panel de CONUS+ es una figura NUEVA -- 2pchi2.f90
nunca produjo un perfil 1D en amplitud, solo la grilla NSI 2D.

Entradas (datos/):
    generic_Xe_perfil.dat      (chi2_nsi_generic + generic_input_Xe.dat)
    generic_conus_perfil.dat   (chi2_nsi_generic + generic_input_conus.dat)
    generic_Xe_resumen.txt, generic_conus_resumen.txt  (A_best, A_90, chi2_min)

Salida: datos/fig_perfil_cruzado_XeConus.png
"""
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()
import matplotlib.ticker as ticker

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
NARANJA = "#a86a43"


def leer_perfil(path):
    filas = []
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.lstrip().startswith("#") or not ln.strip():
                continue
            ln = re.sub(r"(\d)([-+])(\d{2,3})\b", r"\1E\2\3", ln)
            filas.append([float(x) for x in ln.split()])
    a = np.array(filas)
    return a[:, 0], a[:, 1], a[:, 2]

def leer_resumen(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.lstrip().startswith("#") or "=" not in ln:
                continue
            k, v = ln.split("=", 1)
            try:
                out[k.strip()] = float(v.strip())
            except ValueError:
                pass
    return out

def panel(ax, tag, titulo, xmax_extra):
    A, chi2_sin, chi2_con = leer_perfil(f"{BASE}/generic_{tag}_perfil.dat")
    res = leer_resumen(f"{BASE}/generic_{tag}_resumen.txt")
    chi2_min = chi2_con.min()
    A_best = A[np.argmin(chi2_con)]
    A_90 = res["A_90"]

    dchi2_con = chi2_con - chi2_min
    dchi2_sin = chi2_sin - chi2_sin.min()
    # el motor generico escribe col2 == col3 cuando NO hay nuisance de flujo
    # (RED-100/Xe: sigma_alpha <= 0). En ese caso hay un solo perfil.
    hay_nuisance = np.max(np.abs(chi2_con - chi2_sin)) > 1e-6

    if hay_nuisance:
        ax.plot(A, dchi2_con, color="k", lw=1.7, label=r"$\Delta\chi^2$ (con nuisance de flujo)")
        ax.plot(A, dchi2_sin, color="0.55", lw=1.0, ls="--", label=r"$\Delta\chi^2$ (sin nuisance)")
    else:
        ax.plot(A, dchi2_sin, color="k", lw=1.7, label=r"$\Delta\chi^2(A)$")
    ax.axhline(2.706, color=NARANJA, lw=1.2, ls="--", label=r"90% C.L. ($\Delta\chi^2{=}2{,}706$)")
    ax.axhline(1.000, color="#33546e", lw=0.9, ls=":", label=r"$1\sigma$ ($\Delta\chi^2{=}1$)")
    ax.axvline(A_best, color="k", lw=0.7, ls=":", alpha=0.6)
    ax.axvline(A_90, color=NARANJA, lw=0.7, ls=":", alpha=0.8)
    ax.axvline(1.0, color="0.55", lw=0.7, ls=":", alpha=0.5)

    ymax = 9.0
    ax.annotate(f"$A_{{best}}={A_best:.2f}$", xy=(A_best, 0),
                xytext=(A_best + 0.05*xmax_extra, ymax*0.72), fontsize=8.3,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.annotate(f"$A_{{90\\%}}={A_90:.2f}$", xy=(A_90, 2.706),
                xytext=(A_90 - 0.28*xmax_extra, ymax*0.48), fontsize=8.3, color=NARANJA,
                arrowprops=dict(arrowstyle="->", color=NARANJA, lw=0.7))
    ax.text(1.0, ymax*0.05, "SM\n($A{=}1$)", ha="center", va="bottom", fontsize=7.3, color="0.5")

    ax.set_xlim(max(-0.8, A_best - 1.0), A_90 + 0.4*xmax_extra)
    ax.set_ylim(-0.3, ymax)
    ax.set_xlabel(r"Amplitud $A$")
    ax.set_title(titulo)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=6))
    return A_best, A_90, chi2_min

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.4))
Ab_xe, A90_xe, cmin_xe = panel(ax1, "Xe", "Xe (RED-100)", 1.0)
Ab_co, A90_co, cmin_co = panel(ax2, "conus", "CONUS+ (nuevo)", 1.0)
ax1.set_ylabel(r"$\Delta\chi^2 = \chi^2(A) - \chi^2_{\min}$")
ax1.legend(loc="upper center", handlelength=1.6)
ax2.legend(loc="upper left", handlelength=1.6)
fig.suptitle("Perfil $\\Delta\\chi^2(A)$: mismo binario, dos experimentos", y=1.02, fontsize=12)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_perfil_cruzado_XeConus.png")

print("=" * 78)
print(" PERFIL chi2(A) -- MISMO BINARIO chi2_nsi_generic, dos experimentos")
print("=" * 78)
print(f" Xe    : A_best={Ab_xe:.4f}  A_90={A90_xe:.4f}  chi2_min={cmin_xe:.4f}"
      f"   (chi2.f90 real: A_90=107.32, chi2_min=17.51)")
print(f" CONUS+: A_best={Ab_co:.4f}  A_90={A90_co:.4f}  chi2_min={cmin_co:.4f}"
      f"   (2pchi2.f90: chi2_min~7.35)")
print(f"\n  {BASE}/fig_perfil_cruzado_XeConus.png")
