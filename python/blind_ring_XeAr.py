#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
El anillo ciego de la NSI de sola tasa, para Xe y Ar en el MISMO plano.

Panel A -- plano (eps_ee^uV, eps_ee^dV): la carga diagonal es
    q_ee = eps_u (2Z+N) + eps_d (Z+2N),
asi que la "direccion ciega" (q_ee tal que Q_W + q_ee = 0) es una recta
por el origen con angulo theta(Z,N) = arctan[(Z+2N)/(2Z+N)]. Se dibujan
las rectas ciegas de Xe, Ge y Ar; el angulo Xe-Ar (1.44 deg) y el maximo
posible O-Xe (3.35 deg, O-16 con Z=N) se anotan. Es la degeneracion
ESTRUCTURAL: ningun par de blancos reales la rompe apreciablemente.

Panel B -- plano (eps_ee^dV, eps_emu^dV) (ipar=5): el lugar geometrico
q_eff^2 = Q_W^2 es una circunferencia por el SM, centro
(rho_eps, 0) con rho_eps = |Q_W|/(Z+2N) y radio rho_eps. Se superpone el
contorno numerico Delta_chi2 = 4.605 (90%, 2 gdl) de la grilla ideal
(chi2_nsi_2D{Xe,Ar}_ideal.dat): la circunferencia analitica cae
exactamente sobre el borde de la region numerica.

Entradas: datos/chi2_nsi_2D{Xe,Ar}_ideal.dat
Salida  : datos/fig_blind_ring_XeAr.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
S2W = 0.23857
C_XE, C_AR, C_GE = "#33546e", "#a86a43", "#5c7053"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.fontsize": 8.3,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})

TARG = {
    "Xe": dict(Z=54, N=77, c=C_XE),
    "Ge": dict(Z=32, N=41, c=C_GE),
    "Ar": dict(Z=18, N=22, c=C_AR),
}
for name, d in TARG.items():
    Z, N = d["Z"], d["N"]
    d["QW"] = -N / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * Z
    d["ZN2"] = Z + 2 * N
    d["theta"] = np.degrees(np.arctan2(Z + 2 * N, 2 * Z + N))
    d["rho"] = abs(d["QW"]) / d["ZN2"]

th_O = np.degrees(np.arctan2(24, 24))   # O-16: Z=N=8 -> arctan(24/24) = 45 deg

fig, (aA, aB) = plt.subplots(1, 2, figsize=(11.0, 5.2))

# ================================================= Panel A: direccion ciega
t = np.linspace(-1.2, 1.2, 10)
for name, d in TARG.items():
    ang = np.radians(d["theta"])
    aA.plot(t * np.cos(ang), t * np.sin(ang), color=d["c"], lw=2.0,
            label=fr"{name}: $\theta={d['theta']:.2f}^\circ$")
# referencia O-16 (maximo posible)
angO = np.radians(th_O)
aA.plot(t * np.cos(angO), t * np.sin(angO), color="0.55", lw=1.0, ls="--",
        label=fr"$^{{16}}$O (máx.): $\theta={th_O:.2f}^\circ$")
aA.plot(0, 0, "k+", ms=11, mew=1.6)
aA.set_xlim(-1.2, 1.2)
aA.set_ylim(-1.2, 1.2)
aA.set_aspect("equal")
aA.set_xlabel(r"$\varepsilon_{ee}^{uV}$")
aA.set_ylabel(r"$\varepsilon_{ee}^{dV}$")
aA.set_title("(A) Dirección ciega: la degeneración es estructural")
aA.legend(loc="lower right")
aA.annotate(fr"Xe$-$Ar: ${TARG['Xe']['theta']-TARG['Ar']['theta']:.2f}^\circ$",
            xy=(0.62, 0.66), xytext=(0.40, 0.98),
            fontsize=9, ha="left",
            arrowprops=dict(arrowstyle="-", color="0.4", lw=0.7))

# ================================================= Panel B: anillo ciego (ipar=5)
def contour_from_grid(path, level=4.605):
    d = np.loadtxt(path, comments="#")
    x = np.unique(d[:, 0])
    y = np.unique(d[:, 1])
    Z = d[:, 2].reshape(len(y), len(x))
    return x, y, Z

for name in ("Xe", "Ar"):
    d = TARG[name]
    # circunferencia analitica
    ph = np.linspace(0, 2 * np.pi, 400)
    aB.plot(d["rho"] + d["rho"] * np.cos(ph), d["rho"] * np.sin(ph),
            color=d["c"], lw=2.0, ls="--",
            label=fr"{name} analítico ($\rho_\varepsilon={d['rho']:.3f}$)")
    # contorno numerico
    try:
        x, y, Zg = contour_from_grid(f"{BASE}/chi2_nsi_2D{name}_ideal.dat")
        cs = aB.contour(x, y, Zg, levels=[4.605], colors=[d["c"]], linewidths=1.0)
    except Exception as e:
        print(f"  (sin grilla numerica para {name}: {e})")

aB.plot(0, 0, "k+", ms=11, mew=1.6, label="SM")
aB.set_xlim(-0.15, 0.55)
aB.set_ylim(-0.35, 0.35)
aB.set_aspect("equal")
aB.set_xlabel(r"$\varepsilon_{ee}^{dV}$")
aB.set_ylabel(r"$\varepsilon_{e\mu}^{dV}$")
aB.set_title("(B) Anillo ciego: analítico (--) sobre el numérico")
aB.legend(loc="upper right")

fig.suptitle("La NSI de sola tasa es degenerada: un anillo, y ningún par de "
             "blancos rompe la dirección", y=1.00, fontsize=11.5)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_blind_ring_XeAr.png")

print("=" * 68)
for name, d in TARG.items():
    print(f"  {name}: Q_W={d['QW']:8.3f}  theta={d['theta']:.3f} deg  rho_eps={d['rho']:.4f}")
print(f"  Xe-Ar = {TARG['Xe']['theta']-TARG['Ar']['theta']:.2f} deg   "
      f"O-Xe = {th_O-TARG['Xe']['theta']:.2f} deg (max)")
print(f"\n  {BASE}/fig_blind_ring_XeAr.png")
