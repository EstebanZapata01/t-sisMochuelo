#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dsigma/dT del Xe a E_nu=5 MeV con NSI vectorial: la NSI entra solo por
q_eff^2/Q_W^2 = A_amp, una constante independiente de T, asi que
dsigma/dT(eps) = A_amp(eps) * dsigma/dT_SM(T) -- un simple reescalado, sin
distorsion de forma. Cuatro puntos (A_amp verificado con nsi_core):
  eps_ee^dV=-0.05 (A=1.635, diagonal aumentada)
  eps_ee^dV=+0.05 (A=0.520, diagonal reducida)
  eps_emu^dV=+0.10 (A=1.311, cambio de sabor: solo suma, q_ee=0)
  eps_ee^dV=2 eps_0=0.3589 (A=1 EXACTO: la solucion de cancelacion,
    discontinuo sobre el SM -- misma prediccion, eps distinto)

Panel principal: dsigma/dT log-y, con linea vertical punteada en T_max.
Subpanel: cociente con el SM, recortado en T_max (mas alla, las dos
secciones eficaces son cero y el cociente no esta definido).

Entrada: datos/dsigma_bare_Xe.dat (F^2=1, E_nu=5 MeV; dsigma_bare_XeAr.py)
Salida : datos/fig_dsigma_nsi.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, NEGRO, CICLO, coma, formatear_coma
from nsi_core import TARGET, A_amp_grid
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0, "lines.linewidth": 1.3,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
Z, N, QW = TARGET["Xe"]["Z"], TARGET["Xe"]["N"], TARGET["Xe"]["QW"]
eps0 = -QW / (Z + 2 * N)
eps_canc = 2 * eps0

PUNTOS = [
    (rf"$\varepsilon_{{ee}}^{{dV}}={{-}}0{{,}}05$", -0.05, 0.0, "-", CICLO[0], ""),
    (rf"$\varepsilon_{{ee}}^{{dV}}={{+}}0{{,}}05$", 0.05, 0.0, "-", CICLO[1], ""),
    (rf"$\varepsilon_{{e\mu}}^{{dV}}={{+}}0{{,}}10$", 0.0, 0.10, "-", CICLO[2], " (cambio de sabor)"),
    (rf"$\varepsilon_{{ee}}^{{dV}}=2\varepsilon_0={coma(eps_canc,4)}$", eps_canc, 0.0, "--", NEGRO, " (cancelación)"),
]

d = np.loadtxt(f"{BASE}/dsigma_bare_Xe.dat", comments="#")
sel = d[:, 0] == 5.0
T, dsF1 = d[sel, 1], d[sel, 2]
Tmax_keV = d[sel, 4][0]
order = np.argsort(T)
T, dsF1 = T[order], dsF1[order]

fig, (ax, axr) = plt.subplots(2, 1, figsize=(4.60, 4.3),
                                gridspec_kw=dict(height_ratios=[2.2, 1.0], hspace=0.12))
ax.plot(T, dsF1, color=NEGRO, lw=1.8, label="SM ($A=1$)")
axr.plot([0, Tmax_keV], [1.0, 1.0], color=NEGRO, lw=1.4)

print("  A_amp verificados (nsi_core):")
for lab, ex, ey, ls, col, suf in PUNTOS:
    A = A_amp_grid(ex, ey, Z, N, QW, "dV", "dV", "ee", "emu")
    print(f"    {lab}{suf}: A={A:.4f}")
    ax.plot(T, A * dsF1, color=col, ls=ls, lw=1.5, label=lab + rf", $A={coma(A,3)}$" + suf)
    axr.plot([0, Tmax_keV], [A, A], color=col, ls=ls, lw=1.4)

ax.axvline(Tmax_keV, color="0.4", lw=1.0, ls=":")
ax.text(Tmax_keV + 0.012, 3.0e-41, rf"$T^{{\mathrm{{máx}}}}={coma(Tmax_keV,3)}$ keV",
        rotation=90, va="bottom", ha="center", fontsize=8.0, color="0.35")
axr.axvline(Tmax_keV, color="0.4", lw=1.0, ls=":")

ax.set_xlim(0, 0.45)
ax.set_yscale("log")
ax.set_ylabel(r"$d\sigma/dT$  [cm$^2$/keV]")
ax.legend(loc="lower left", fontsize=8.0, framealpha=0.9)
plt.setp(ax.get_xticklabels(), visible=False)

axr.set_xlim(0, 0.45)
axr.set_ylabel(r"$(d\sigma/dT)_{\rm NSI}\,/\,(d\sigma/dT)_{\rm SM}$", fontsize=8.5)
axr.set_xlabel(r"Energía de retroceso nuclear, $T$  [keV]")

formatear_coma(ax, x=False, y=False)
formatear_coma(axr)

fig.tight_layout()
out = f"{BASE}/fig_dsigma_nsi.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
dest = f"{DEST}/fig_dsigma_nsi.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"\n  T_max(5 MeV, Xe) = {Tmax_keV:.4f} keV")
print(f"  -> {dest}  (md5 {md5(out)})")
