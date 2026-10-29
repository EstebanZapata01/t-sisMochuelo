#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONUS+ (Phys. Rev. D 111, 075025), 4 planos (chi2_nsi_2D, ipar 1/2/5/7),
regiones 68/90/95% (Dchi2=2,30/4,61/5,99 sobre el chi2_min global) rellenas
en tonos anidados de C_GE (68% el mas oscuro, 95% el mas claro) con los
contornos punteado/solido/raya-punto encima. En trazo fino, la curva
analitica A_amp(eps)=A+-(90%) (de generic_conus_perfil.dat, perfilando
alpha): como chi2(eps) solo depende de A_amp=q_eff^2/Q_W^2, coincide EXACTO
con el contorno numerico del 90% (se verifica en stdout). En (c) y (d) se
marca tambien el centro (eps_0,0), A_amp=0 (tasa nula, excluido).

Dibujada ya al tamano final de pagina (ancho completo, 6,14 in).

Entrada: datos/chi2_nsi_2Dconus{,_ipar1,_ipar2,_ipar7}.dat,
         datos/generic_conus_perfil.dat
Salida : datos/fig_conus_planos.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_GE, NEGRO, GRIS, coma, formatear_coma, nota
from nsi_core import A_amp_grid, COEF
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.5, "lines.linewidth": 1.15,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"

Z_Ge, N_Ge = 32.0, 41.0
S2W = 0.23857
QW_SM = -N_Ge / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * Z_Ge
eps0_Ge = -QW_SM / COEF["dV"](Z_Ge, N_Ge)

d = np.loadtxt(f"{BASE}/generic_conus_perfil.dat", comments="#")
Ap, chi2p = d[:, 0], d[:, 2]
chi2min = chi2p.min()
dchi2p = chi2p - chi2min


def A_roots(level):
    f = lambda A: np.interp(A, Ap, dchi2p) - level
    i0 = np.argmin(dchi2p)
    Alo = brentq(f, Ap[0] + 1e-6, Ap[i0])
    Ahi = brentq(f, Ap[i0], Ap[-1] - 1e-6)
    return Alo, Ahi


NIVELES = [(2.30, "68%", ":"), (4.61, "90%", "-"), (5.99, "95%", "-.")]
ROOTS = {lvl: A_roots(lvl) for lvl, _, _ in NIVELES}
A90m, A90p = ROOTS[4.61]
ALPHAS = [0.55, 0.32, 0.16]   # 68% (adentro, mas oscuro) -> 95% (afuera, mas claro)

PANELES = [
    (1, "dV", "uV", "ee", "ee", r"$\varepsilon_{ee}^{dV}$", r"$\varepsilon_{ee}^{uV}$",
     (-1, 1), (-1, 1), "(a) dos bandas", False),
    (2, "dV", "uV", "emu", "emu", r"$\varepsilon_{e\mu}^{dV}$", r"$\varepsilon_{e\mu}^{uV}$",
     (-1, 1), (-1, 1), "(b) una banda", False),
    (5, "dV", "dV", "ee", "emu", r"$\varepsilon_{ee}^{dV}$", r"$\varepsilon_{e\mu}^{dV}$",
     (-0.12, 0.48), (-0.3, 0.3), "(c) anillo", True),
    (7, "dV", "uV", "ee", "emu", r"$\varepsilon_{ee}^{dV}$", r"$\varepsilon_{e\mu}^{uV}$",
     (-0.12, 0.48), (-0.3, 0.3), "(d) anillo elíptico", True),
]

fig, axs = plt.subplots(2, 2, figsize=(6.14, 6.1))
for ax, (ipar, fx, fy, lx, ly, xl, yl, xlim, ylim, titulo, marca_centro) in zip(axs.flat, PANELES):
    fname = f"{BASE}/chi2_nsi_2Dconus.dat" if ipar == 5 else f"{BASE}/chi2_nsi_2Dconus_ipar{ipar}.dat"
    dd = np.loadtxt(fname)
    n = int(round(np.sqrt(dd.shape[0])))
    X = dd[:, 0].reshape(n, n); Y = dd[:, 1].reshape(n, n); Zc = dd[:, 2].reshape(n, n)
    dchi2 = Zc - chi2min

    niveles_fill = [0.0] + [lvl for lvl, _, _ in NIVELES]
    colores_fill = [(*mpl.colors.to_rgb(C_GE), a) for a in ALPHAS]
    ax.contourf(X, Y, dchi2, levels=niveles_fill, colors=colores_fill)
    for lvl, lab, ls in NIVELES:
        ax.contour(X, Y, dchi2, levels=[lvl], colors=[C_GE], linewidths=1.0, linestyles=[ls])

    g = np.linspace(xlim[0], xlim[1], 500)
    EX, EY = np.meshgrid(g, np.linspace(ylim[0], ylim[1], 500))
    Amap = A_amp_grid(EX, EY, Z_Ge, N_Ge, QW_SM, fx, fy, lx, ly)
    ax.contour(EX, EY, Amap, levels=[A90m, A90p], colors=["0.15"], linewidths=0.6, linestyles=["-"])
    ax.contour(EX, EY, Amap, levels=[1.0], colors=[GRIS], linewidths=0.7, linestyles=[":"])

    ax.plot(0, 0, "*", ms=8, color=NEGRO, zorder=6)
    if marca_centro:
        ax.plot(eps0_Ge, 0, "x", ms=6, mew=1.4, color=NEGRO, zorder=6)
        ax.text(eps0_Ge, -0.035, r"$A_{\rm amp}=0$", ha="center", va="top",
                color=NEGRO, fontsize=7.2)
    ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_aspect("equal")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_title(titulo, fontsize=9.5)
    formatear_coma(ax)

axs[0, 0].plot([], [], color=C_GE, lw=1.0, ls=":", label="68%")
axs[0, 0].plot([], [], color=C_GE, lw=1.0, ls="-", label="90%")
axs[0, 0].plot([], [], color=C_GE, lw=1.0, ls="-.", label="95%")
axs[0, 0].plot([], [], color="0.15", lw=0.6, label=r"$A_{\rm amp}=A_\pm$ (90%, analítico)")
axs[0, 0].plot([], [], color=GRIS, lw=0.7, ls=":", label=r"$A_{\rm amp}=1$ (SM)")
axs[0, 0].legend(loc="lower left", fontsize=6.8)

fig.tight_layout()
out = f"{BASE}/fig_conus_planos.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = f"{DEST}/fig_conus_planos.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print("  Verificacion contorno numerico vs curva analitica A_amp=A+- (90%):")
print("  (por construccion exacta: chi2(eps) depende solo de A_amp(eps), ver nsi_core)")
for lvl, lab, _ in NIVELES:
    Alo, Ahi = ROOTS[lvl]
    print(f"    Dchi2={lvl} ({lab}): A-={Alo:.4f}  A+={Ahi:.4f}")
print(f"  eps_0(Ge) = {eps0_Ge:.4f}")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
