#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Factor de forma nuclear de Helm F(q^2) para Xe y Ar, y cuanto se aparta
de la aproximacion F^2 = 1 que usa el codigo (xsections_nest.f90).

    F(q^2) = 3 j_1(qR_A)/(qR_A) * exp(-q^2 s^2 / 2)
    R_A ~ 1.2 A^{1/3} fm ,  s ~ 0.9 fm ,  q = sqrt(2 M T)

El documento (Ec. helm) afirma EN TEXTO que qR_A ~ 0.1 y F^2 >= 0.99 en la
ROI. Eso NO es correcto: para Xe a T = 1 keV, qR_A ~ 0.48 y F^2 ~ 0.95.
Este script cuantifica el tamano REAL de la aproximacion F^2 = 1: ~5% en
la tasa absoluta de Xe y ~4% en el cociente F^2_Xe/F^2_Ar. La direccion
es conservadora: incluir Helm SUPRIME mas al Xe que al Ar, asi que la
ventaja de Ar y el A_90 de Xe (107) se moverian, si acaso, a favor de
las conclusiones ya sacadas. La geometria NSI (anillo ciego, angulo
1.44 deg) no cambia: F^2 es un factor global, no afecta la direccion
ciega (que depende de Q_W).

Salida: datos/fig_helm_FF.png  (+ numeros a stdout)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
HBARC = 197.3269804        # MeV*fm
AMU   = 931.49410242       # MeV
S_FM  = 0.9                # fm
C_XE, C_AR = "#33546e", "#a86a43"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.fontsize": 8.5,
    "lines.linewidth": 1.8,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})


def helm_F2(T_keV, A):
    """F^2 de Helm vs energia de retroceso T [keV] para masa A."""
    M = A * AMU                              # MeV
    T = T_keV * 1e-3                         # MeV
    q = np.sqrt(2.0 * M * T)                 # MeV
    R = 1.2 * A ** (1.0 / 3.0)              # fm
    x = q * R / HBARC                        # adimensional
    qs = q * S_FM / HBARC
    j1 = np.sin(x) / x**2 - np.cos(x) / x
    F = 3.0 * j1 / x * np.exp(-0.5 * qs**2)
    return F**2


T = np.linspace(0.05, 6.0, 2000)
F2_xe = helm_F2(T, 131)
F2_ar = helm_F2(T, 40)

# ROI CEvNS (retroceso nuclear) aproximada para cada blanco
# Xe: N_e = 4-7 -> T ~ 0.35-0.65 keV (metodologia.tex Sec. cinematica /
# threshold_table.py: T_thr(N_e>=4) ~ 0.64 keV). Ar: N_e = 1-5.
ROI_XE = (0.35, 0.65)
ROI_AR = (0.15, 3.5)

def in_roi(Tarr, F2, roi):
    m = (Tarr >= roi[0]) & (Tarr <= roi[1])
    return F2[m].min(), F2[m].max()

xe_lo, xe_hi = in_roi(T, F2_xe, ROI_XE)
ar_lo, ar_hi = in_roi(T, F2_ar, ROI_AR)

print("=" * 70)
print(" Factor de forma de Helm en la ROI CEvNS")
print("=" * 70)
print(f"  Xe (T in {ROI_XE} keV):  F^2 = {xe_lo:.4f} .. {xe_hi:.4f}")
print(f"  Ar (T in {ROI_AR} keV):  F^2 = {ar_lo:.4f} .. {ar_hi:.4f}")
_troi = np.linspace(ROI_XE[0], ROI_XE[1], 50)
ratio_roi = helm_F2(_troi, 131) / helm_F2(_troi, 40)
print(f"  cociente F^2_Xe/F^2_Ar en la ROI de Xe: {ratio_roi.min():.4f} .. {ratio_roi.max():.4f}")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.2, 4.2))

a1.plot(T, F2_xe, color=C_XE, label="Xe ($A=131$)")
a1.plot(T, F2_ar, color=C_AR, label="Ar ($A=40$)")
a1.axhline(1.0, color="0.5", lw=0.8, ls="--", label=r"aproximación $F^2=1$")
a1.axvspan(*ROI_XE, color=C_XE, alpha=0.10)
a1.axvspan(*ROI_AR, color=C_AR, alpha=0.08)
a1.set_xlabel(r"$T_{\rm nr}$ [keV]")
a1.set_ylabel(r"$F^2(q^2)$")
a1.set_title("Factor de forma de Helm vs. retroceso")
a1.set_xlim(0, 6)
a1.set_ylim(0.90, 1.005)
a1.legend(loc="lower left")
a1.text(0.97, 0.05, "bandas: ROI CEνNS\nde cada blanco",
        transform=a1.transAxes, ha="right", va="bottom", fontsize=7.8, color="0.4")

r = F2_xe / F2_ar
a2.plot(T, r, color="black")
a2.axhline(1.0, color="0.5", lw=0.8, ls="--")
a2.axvspan(*ROI_XE, color="0.88")
a2.set_xlabel(r"$T_{\rm nr}$ [keV]")
a2.set_ylabel(r"$F^2_{\rm Xe}\,/\,F^2_{\rm Ar}$")
a2.set_title(r"El cociente $\ne 1$: la aproximación $F^2=1$ favorece a Xe")
a2.set_xlim(0, 6)
a2.set_ylim(0.90, 1.02)
a2.text(0.5, 0.06, "ROI Xe", transform=a2.get_xaxis_transform(),
        ha="center", fontsize=7.8, color="0.35")
a2.text(0.97, 0.10,
        rf"en la ROI de Xe: $F^2_{{\rm Xe}}/F^2_{{\rm Ar}}\approx{ratio_roi.mean():.2f}$"
        "\n(incluir Helm refuerza la ventaja de Ar)",
        transform=a2.transAxes, ha="right", va="bottom", fontsize=7.8)

fig.suptitle("La aproximación $F^2=1$ cuesta ~5% en la tasa de Xe y ~4% en el "
             "cociente Xe/Ar, en dirección conservadora", y=1.02, fontsize=11)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_helm_FF.png")
print(f"\n  {BASE}/fig_helm_FF.png")
