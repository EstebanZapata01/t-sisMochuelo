#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
El anillo ciego de la NSI de sola tasa, para Xe y Ar en el MISMO plano.

Panel A -- caracter ESTRUCTURAL de la direccion ciega. El angulo ciego es
    theta = arctan f(r),   f(r) = (1+2r)/(2+r),   r = N/Z,
monotona creciente y acotada en [45 deg, arctan(2)=63.43 deg]. Todos los
nucleos estables viven en el tramo casi plano (r ~ 1-1.6), de modo que
ningun par de blancos CE$\nu$NS reales separa sus direcciones ciegas mas
de ~3 deg. Inset: el zoom con Delta_phi(Xe,Ar) = 1.44 deg.

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
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO, GRIS, NEGRO, zoom, nota, sombra
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
S2W = 0.23857

TARG = {
    "$^{16}$O":  dict(Z=8,  N=8.0,  c="0.55"),
    "$^{23}$Na": dict(Z=11, N=12.0, c="0.55"),
    "$^{40}$Ar": dict(Z=18, N=22.0, c=C_AR),
    "$^{73}$Ge": dict(Z=32, N=41.0, c=C_GE),
    "$^{127}$I": dict(Z=53, N=74.0, c="0.55"),
    "$^{131}$Xe":dict(Z=54, N=131.293 - 54, c=C_XE),  # N=<A>-Z, mezcla natural
}
for name, d in TARG.items():
    Z, N = d["Z"], d["N"]
    d["QW"] = -N / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * Z
    d["ZN2"] = Z + 2 * N
    d["r"] = N / Z
    d["theta"] = np.degrees(np.arctan2(Z + 2 * N, 2 * Z + N))
    d["rho"] = abs(d["QW"]) / d["ZN2"]

RING = {"Xe": TARG["$^{131}$Xe"], "Ar": TARG["$^{40}$Ar"]}
th_xe, th_ar = TARG["$^{131}$Xe"]["theta"], TARG["$^{40}$Ar"]["theta"]
th_O = TARG["$^{16}$O"]["theta"]

fig, (aA, aB) = plt.subplots(1, 2, figsize=(11.2, 5.0))

def theta_f(rr):
    return np.degrees(np.arctan((1.0 + 2.0 * rr) / (2.0 + rr)))

# ============================== Panel A: theta = arctan f(r), estructural
r = np.linspace(0.9, 3.2, 500)
aA.plot(r, theta_f(r), color=NEGRO, lw=1.8)
aA.axhline(theta_f(1e6), color=GRIS, ls="--", lw=0.8)
aA.text(3.15, theta_f(1e6) - 1.4, r"$\theta\to 63{,}4^\circ$", ha="right",
        fontsize=8, color=GRIS)
sombra(aA, 1.0, 1.6, color=C_XE, alpha=0.09)
aA.text(1.30, 41.5, "valle de\nestabilidad", ha="center", fontsize=8, color=C_XE)
for name, d in TARG.items():
    aA.plot(d["r"], d["theta"], "o", ms=4.5, color=d["c"], zorder=5)
aA.set_xlim(0.9, 3.2)
aA.set_ylim(40, 65)
aA.set_xlabel(r"$r = N/Z$")
aA.set_ylabel(r"$\theta(r) = \arctan\dfrac{1+2r}{2+r}$   [$^\circ$]")

# --- inset: zoom al cúmulo de blancos reales, con Delta_phi Xe-Ar ---
axA = zoom(aA, [0.50, 0.08, 0.46, 0.46], (1.16, 1.52), (46.6, 48.9))
axA.plot(r, theta_f(r), color=NEGRO, lw=1.4)
for key, col in (("$^{40}$Ar", C_AR), ("$^{73}$Ge", C_GE), ("$^{131}$Xe", C_XE)):
    d = TARG[key]
    axA.plot(d["r"], d["theta"], "o", ms=5, color=col)
    axA.axhline(d["theta"], color=col, lw=0.5, ls=":")
axA.annotate("", xy=(1.185, th_xe), xytext=(1.185, th_ar),
             arrowprops=dict(arrowstyle="<->", lw=1.0, color=NEGRO))
axA.text(1.20, (th_xe + th_ar) / 2,
         rf"$\Delta\phi_{{\rm Xe,Ar}} = {th_xe - th_ar:.2f}^\circ$",
         va="center", fontsize=8.2)

# ================================================= Panel B: anillo ciego (ipar=5)
def contour_from_grid(path):
    d = np.loadtxt(path, comments="#")
    x = np.unique(d[:, 0]); y = np.unique(d[:, 1])
    return x, y, d[:, 2].reshape(len(y), len(x))

# Delta_chi2 critico (2 g.d.l.) para 68%, 90% y 95% C.L.
CL_LEVELS = [(2.30, "68%", ":"), (4.605, "90%", "-"), (5.99, "95%", "-.")]

ph = np.linspace(0, 2 * np.pi, 512)
def dibuja_anillo(ax, name, lw_a=1.8):
    d = RING[name]
    ax.plot(d["rho"] * (1 + np.cos(ph)), d["rho"] * np.sin(ph),
            color=d["c"], lw=lw_a, ls=(0, (5, 2)))
    try:
        x, y, Zg = contour_from_grid(f"{BASE}/chi2_nsi_2D{name}_ideal.dat")
        for dchi2, _, ls in CL_LEVELS:
            ax.contour(x, y, Zg, levels=[dchi2], colors=[d["c"]],
                       linewidths=0.9, linestyles=[ls])
    except Exception as e:
        print(f"  (sin grilla numerica para {name}: {e})")

for name in ("Xe", "Ar"):
    dibuja_anillo(aB, name)
    aB.plot([], [], color=RING[name]["c"], lw=1.8, ls=(0, (5, 2)),
            label=fr"{name}: $\rho_\varepsilon = {RING[name]['rho']:.3f}$")
for _, lab, ls in CL_LEVELS:
    aB.plot([], [], color=GRIS, lw=0.9, ls=ls, label=lab + " C.L.")
aB.plot(0, 0, "+", ms=10, mew=1.4, color=NEGRO)
aB.set_xlim(-0.12, 0.52)
aB.set_ylim(-0.32, 0.32)
aB.set_aspect("equal")
aB.set_xlabel(r"$\varepsilon_{ee}^{dV}$")
aB.set_ylabel(r"$\varepsilon_{e\mu}^{dV}$")
aB.legend(loc="upper left", ncol=2, fontsize=7.6)

# --- inset: el borde del anillo, analítico (--) sobre el numérico ---
axB = zoom(aB, [0.60, 0.06, 0.38, 0.38], (0.30, 0.40), (-0.05, 0.05))
for name in ("Xe", "Ar"):
    dibuja_anillo(axB, name, lw_a=1.6)
axB.set_aspect("equal")

fig.tight_layout()
fig.savefig(f"{BASE}/fig_blind_ring_XeAr.png")

dphi = th_xe - th_ar
print("=" * 68)
for name, d in TARG.items():
    print(f"  {name:10s}: r={d['r']:.3f}  theta={d['theta']:.3f} deg  rho_eps={d['rho']:.4f}")
print(f"  Delta_phi(Xe,Ar) = {dphi:.3f} deg   1/sin = {1/np.sin(np.radians(dphi)):.1f}")
print(f"  O-Xe = {th_xe-th_O:.2f} deg (max)   1/sin = {1/np.sin(np.radians(th_xe-th_O)):.1f}")
print(f"\n  {BASE}/fig_blind_ring_XeAr.png")
