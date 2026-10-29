#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lugares ciegos de Xe y Ar en 4 familias de planos, comparacion IDEAL
(Asimov, 192 kg dia, N_e>=1, F de NEST, sin fondo). Reemplaza
fig_blind_rings.png, fig_lugares_ciegos_familias.png y fig_blind_theta.png.

(a) rectas (eps_ee^dV, eps_ee^uV) y (c) elipses (eps_ee^dV, eps_emu^uV):
    bandas 68/90/95% (Dchi2=2,30/4,61/5,99) de Xe y Ar, analiticas --
    chi2(eps) = (1-A_amp(eps))^2 N_tot es cuadratica en A_amp, asi que
    los bordes son A_amp = 1 +- sqrt(Dchi2/N_tot), cerrados.
(b) circunferencias (eps_ee^dV, eps_emu^dV): igual que fig_blind_rings
    (O, Na, Ge, I en gris; Xe/Ar con el contorno NUMERICO de
    chi2_nsi_2D{Xe,Ar}_ideal.dat, Dchi2=2,30/4,61/5,99, mas el anillo
    analitico Dchi2=0 de referencia).
(d) theta = arctan[(1+2r)/(2+r)] vs r=N/Z (igual que fig_blind_theta),
    con el recuadro Xe-Ar.

Entrada: datos/chi2_nsi_2D{Xe,Ar}_ideal.dat, datos/sensib_ideal_{Xe,Ar}.dat
Salida : datos/fig_lugares_ciegos_XeAr.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, NEGRO, GRIS, zoom, coma, formatear_coma, nota
from nsi_core import TARGET, A_amp_grid, COEF
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.5, "lines.linewidth": 1.15,
})

def halo_ticks(axin):
    """Caja blanca detras de los numeros del recuadro (si una curva pasa justo
    detras, si no quedan tapados)."""
    for lbl in list(axin.get_xticklabels()) + list(axin.get_yticklabels()):
        lbl.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.5))


BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
S2W = 0.23857
COL = {"Xe": C_XE, "Ar": C_AR}
NIVELES = [(2.30, "68%", ":"), (4.61, "90%", "-"), (5.99, "95%", "-.")]

# N_tot = exposicion(192) * R_tot(N_e>=1), F de NEST (sensib_ideal_{Xe,Ar}.dat)
N_tot = {}
for tag in ("Xe", "Ar"):
    d = np.loadtxt(f"{BASE}/sensib_ideal_{tag}.dat", comments="#", dtype=str)
    fila = d[(d[:, 0] == "1") & (d[:, 1] == "Fnest")][0]
    N_tot[tag] = float(fila[4])   # Ntot_ROI@192

fig, axs = plt.subplots(2, 2, figsize=(6.14, 6.1))

# ======================================================= (a) y (c): bandas analiticas
def panel_bandas(ax, fx, fy, lx, ly, xlim, ylim, titulo):
    g = np.linspace(*xlim, 700)
    EX, EY = np.meshgrid(g, np.linspace(*ylim, 700))
    for tag in ("Xe", "Ar"):
        t = TARGET[tag]
        Amap = A_amp_grid(EX, EY, t["Z"], t["N"], t["QW"], fx, fy, lx, ly)
        for lvl, lab, ls in NIVELES:
            Am, Ap = 1.0 - np.sqrt(lvl / N_tot[tag]), 1.0 + np.sqrt(lvl / N_tot[tag])
            ax.contour(EX, EY, Amap, levels=[Am, Ap], colors=[COL[tag]], linewidths=1.2, linestyles=[ls])
        ax.plot([], [], color=COL[tag], lw=1.2, label=tag)
    ax.plot(0, 0, "*", ms=10, color=NEGRO, zorder=6)
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_title(titulo, fontsize=10)
    formatear_coma(ax)
    return g


def sin_numeros(axin):
    """Solo las marcas del inserto, sin numeros: el pie ya dice que amplia
    cada recuadro, y los numeros (a la derecha) terminaban saliendose del
    marco del panel padre o cortados por el borde de la figura."""
    axin.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)


panel_bandas(axs[0, 0], "dV", "uV", "ee", "ee", (-1, 1), (-1, 1), "(a) rectas")
axs[0, 0].set_xlabel(r"$\varepsilon_{ee}^{dV}$"); axs[0, 0].set_ylabel(r"$\varepsilon_{ee}^{uV}$")
zc = zoom(axs[0, 0], [0.56, 0.06, 0.40, 0.34], (0.28, 0.42), (-0.06, 0.06))
g2 = np.linspace(0.28, 0.42, 300); EY2g = np.linspace(-0.06, 0.06, 300)
EX2, EY2 = np.meshgrid(g2, EY2g)
for tag in ("Xe", "Ar"):
    t = TARGET[tag]
    Amap = A_amp_grid(EX2, EY2, t["Z"], t["N"], t["QW"], "dV", "uV", "ee", "ee")
    for lvl, lab, ls in NIVELES:
        Am, Ap = 1.0 - np.sqrt(lvl / N_tot[tag]), 1.0 + np.sqrt(lvl / N_tot[tag])
        zc.contour(EX2, EY2, Amap, levels=[Am, Ap], colors=[COL[tag]], linewidths=1.1, linestyles=[ls])
formatear_coma(zc)
sin_numeros(zc)
halo_ticks(zc)

panel_bandas(axs[1, 0], "dV", "uV", "ee", "emu", (-0.05, 0.52), (-0.3, 0.3), "(c) elipses")
axs[1, 0].set_xlabel(r"$\varepsilon_{ee}^{dV}$"); axs[1, 0].set_ylabel(r"$\varepsilon_{e\mu}^{uV}$")
axs[1, 0].set_aspect("equal")
zc2 = zoom(axs[1, 0], [0.56, 0.06, 0.40, 0.34], (0.28, 0.42), (-0.06, 0.06))
EX3, EY3 = np.meshgrid(g2, EY2g)
for tag in ("Xe", "Ar"):
    t = TARGET[tag]
    Amap = A_amp_grid(EX3, EY3, t["Z"], t["N"], t["QW"], "dV", "uV", "ee", "emu")
    for lvl, lab, ls in NIVELES:
        Am, Ap = 1.0 - np.sqrt(lvl / N_tot[tag]), 1.0 + np.sqrt(lvl / N_tot[tag])
        zc2.contour(EX3, EY3, Amap, levels=[Am, Ap], colors=[COL[tag]], linewidths=1.1, linestyles=[ls])
formatear_coma(zc2)
sin_numeros(zc2)
halo_ticks(zc2)

# ======================================================= (b) circunferencias (= fig_blind_rings)
ax = axs[0, 1]
TARG = {
    "$^{16}$O":  dict(Z=8,  N=8.0,  c="0.55"),
    "$^{23}$Na": dict(Z=11, N=12.0, c="0.55"),
    "$^{40}$Ar": dict(Z=18, N=22.0, c=C_AR),
    "$^{73}$Ge": dict(Z=32, N=41.0, c=C_GE),
    "$^{127}$I": dict(Z=53, N=74.0, c="0.55"),
    "Xe nat.":   dict(Z=54, N=131.3879 - 54, c=C_XE),
}
for name, d in TARG.items():
    Z, N = d["Z"], d["N"]
    d["QW"] = -N / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * Z
    d["rho"] = abs(d["QW"]) / (Z + 2 * N)
RING = {"Xe": TARG["Xe nat."], "Ar": TARG["$^{40}$Ar"]}
ph = np.linspace(0, 2 * np.pi, 512)
ESTRUCT = ("$^{16}$O", "$^{23}$Na", "$^{73}$Ge", "$^{127}$I")
for name in ESTRUCT:
    d = TARG[name]
    ax.plot(d["rho"] * (1 + np.cos(ph)), d["rho"] * np.sin(ph), color=GRIS, lw=0.7, alpha=0.8, zorder=1)


def contour_grid(path):
    d = np.loadtxt(path, comments="#")
    x = np.unique(d[:, 0]); y = np.unique(d[:, 1])
    return x, y, d[:, 2].reshape(len(y), len(x))


for tag in ("Xe", "Ar"):
    d = RING[tag]
    ax.plot(d["rho"] * (1 + np.cos(ph)), d["rho"] * np.sin(ph), color=COL[tag], lw=1.6, ls=(0, (5, 2)))
    x, y, Zg = contour_grid(f"{BASE}/chi2_nsi_2D{tag}_ideal.dat")
    for lvl, lab, ls in NIVELES:
        ax.contour(x, y, Zg, levels=[lvl], colors=[COL[tag]], linewidths=1.0, linestyles=[ls])
ax.plot(0, 0, "*", ms=8, color=NEGRO, zorder=6)
ax.set_xlim(-0.06, 0.52); ax.set_ylim(-0.3, 0.3)
ax.set_aspect("equal")
ax.set_xlabel(r"$\varepsilon_{ee}^{dV}$"); ax.set_ylabel(r"$\varepsilon_{e\mu}^{dV}$")
ax.set_title("(b) circunferencias", fontsize=10)
formatear_coma(ax)
axz = zoom(ax, [0.56, 0.06, 0.38, 0.34], (0.30, 0.40), (-0.05, 0.05))
for name in ESTRUCT:
    d = TARG[name]
    axz.plot(d["rho"] * (1 + np.cos(ph)), d["rho"] * np.sin(ph), color=GRIS, lw=0.7, alpha=0.8, zorder=1)
for tag in ("Xe", "Ar"):
    d = RING[tag]
    axz.plot(d["rho"] * (1 + np.cos(ph)), d["rho"] * np.sin(ph), color=COL[tag], lw=1.4, ls=(0, (5, 2)))
    x, y, Zg = contour_grid(f"{BASE}/chi2_nsi_2D{tag}_ideal.dat")
    for lvl, lab, ls in NIVELES:
        axz.contour(x, y, Zg, levels=[lvl], colors=[COL[tag]], linewidths=0.9, linestyles=[ls])
axz.set_aspect("equal")
formatear_coma(axz)
sin_numeros(axz)
halo_ticks(axz)

# ======================================================= (d) theta vs r (= fig_blind_theta)
ax = axs[1, 1]


def theta_f(rr):
    return np.degrees(np.arctan((1.0 + 2.0 * rr) / (2.0 + rr)))


r = np.linspace(0.9, 3.2, 500)
ax.plot(r, theta_f(r), color=NEGRO, lw=1.3)
ax.axhline(theta_f(1e6), color=GRIS, ls="--", lw=0.7)
ax.text(3.15, theta_f(1e6) - 1.8, r"$\theta\to 63{,}4^\circ$", ha="right", fontsize=7.2, color=GRIS)
OFF_D = {"$^{16}$O": (5, -11), "$^{23}$Na": (5, -11), "$^{40}$Ar": (-7, 7), "$^{73}$Ge": (5, 5),
         "$^{127}$I": (-10, 11), "Xe nat.": (5, -12)}
for name, d in TARG.items():
    d["r"] = d["N"] / d["Z"]
    d["theta"] = np.degrees(np.arctan2(d["Z"] + 2 * d["N"], 2 * d["Z"] + d["N"]))
    ax.plot(d["r"], d["theta"], "o", ms=5, color=d["c"], zorder=5)
    if name not in ("$^{40}$Ar", "$^{73}$Ge", "Xe nat."):
        ax.annotate(name, (d["r"], d["theta"]), textcoords="offset points",
                    xytext=OFF_D[name], fontsize=7, ha="center")
ax.set_xlim(0.9, 3.2); ax.set_ylim(40, 65)
ax.set_xlabel(r"$r = N/Z$")
ax.set_ylabel(r"$\theta(r)$  [$^\circ$]")
ax.set_title("(d) dirección ciega", fontsize=9.5)
formatear_coma(ax)
axin = zoom(ax, [0.44, 0.06, 0.52, 0.36], (1.16, 1.52), (46.6, 48.9))
axin.plot(r, theta_f(r), color=NEGRO, lw=1.1)
th_xe, th_ar = TARG["Xe nat."]["theta"], TARG["$^{40}$Ar"]["theta"]
IN_X0, IN_X1, IN_Y0, IN_Y1 = 1.16, 1.52, 46.6, 48.9
fx = lambda f: IN_X0 + (IN_X1 - IN_X0) * f
fy = lambda f: IN_Y0 + (IN_Y1 - IN_Y0) * f
for key, col in (("$^{40}$Ar", C_AR), ("$^{73}$Ge", C_GE), ("Xe nat.", C_XE)):
    d = TARG[key]
    axin.plot(d["r"], d["theta"], "o", ms=5, color=col)
axin.text(fx(0.78), fy(0.68), "Xe nat.", ha="left", va="center", fontsize=6.6, color=C_XE)
axin.annotate("$^{73}$Ge", (TARG["$^{73}$Ge"]["r"], TARG["$^{73}$Ge"]["theta"]),
              textcoords="offset points", xytext=(6, -10), ha="left", va="top",
              fontsize=6.6, color=C_GE)
axin.text(fx(0.2), TARG["$^{40}$Ar"]["theta"], "$^{40}$Ar", ha="left", va="center",
          fontsize=6.6, color=C_AR)
axin.annotate("", xy=(1.185, th_xe), xytext=(1.185, th_ar),
              arrowprops=dict(arrowstyle="<->", lw=0.9, color=NEGRO))
axin.text(fx(0.12), fy(0.6), rf"$\Delta\theta={coma(th_xe - th_ar)}^\circ$",
          va="center", ha="left", fontsize=6.8)
formatear_coma(axin)
sin_numeros(axin)
halo_ticks(axin)

# ======================================================= leyenda unica de la figura
from matplotlib.lines import Line2D
leg_handles = [
    Line2D([], [], color=C_XE, lw=1.4, label="Xe"),
    Line2D([], [], color=C_AR, lw=1.4, label="Ar"),
    Line2D([], [], color="0.3", lw=1.1, ls=":", label="68%"),
    Line2D([], [], color="0.3", lw=1.1, ls="-", label="90%"),
    Line2D([], [], color="0.3", lw=1.1, ls="-.", label="95%"),
    Line2D([], [], color=GRIS, lw=0.9, label="otros núcleos (O, Na, Ge, I)"),
]
fig.tight_layout()
leg = fig.legend(handles=leg_handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
                  ncol=3, fontsize=8.0, frameon=False, columnspacing=1.3)
out = f"{BASE}/fig_lugares_ciegos_XeAr.png"
fig.savefig(out, dpi=300, bbox_extra_artists=(leg,), bbox_inches="tight")
plt.close(fig)
dest = f"{DEST}/fig_lugares_ciegos_XeAr.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print(f"  N_tot (192 kgd, N_e>=1, F nest): Xe={N_tot['Xe']:.2f}  Ar={N_tot['Ar']:.2f}")
print(f"  Delta_theta(Xe,Ar) = {th_xe-th_ar:.3f} deg")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
