#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONUS+ (Phys. Rev. D 111, 075025): del perfil en A_amp al intervalo en
epsilon_ee^dV de Ge. Ancho completo de pagina, dos paneles.

(a) Dchi2(A_amp) perfilando alpha (generic_conus_perfil.dat), con horizontales
    en Dchi2=1/2,30/4,61/5,99 y las raices A- , A+ de cada nivel marcadas.
(b) A_amp(eps) = [1+(Z+2N)eps/Q_W]^2 para Ge (Z=32, N=41): la banda
    Dchi2<=1 de (a), trasladada al eje eps por esta parabola, da DOS
    intervalos (rama cerca de eps=0 y rama especular). Se verifica en
    stdout que coinciden con los publicados por CONUS+.

Entrada: datos/generic_conus_perfil.dat
Salida : datos/fig_conus_A_eps.png (+ copia en figuras_metodologia)
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
from nsi_core import COEF
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
QW = -N_Ge / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * Z_Ge
COEF_DV = COEF["dV"](Z_Ge, N_Ge)
PUBLICADO = [(-0.034, 0.024), (0.322, 0.38)]   # CONUS+, PRD 111, 075025 (cifras_tesis_helm.md #5)

d = np.loadtxt(f"{BASE}/generic_conus_perfil.dat", comments="#")
Ap, chi2p = d[:, 0], d[:, 2]
chi2min = chi2p.min()
dchi2p = chi2p - chi2min
i0 = np.argmin(dchi2p)
A_min = Ap[i0]


def A_roots(level):
    f = lambda A: np.interp(A, Ap, dchi2p) - level
    Alo = brentq(f, Ap[0] + 1e-9, Ap[i0])
    Ahi = brentq(f, Ap[i0], Ap[-1] - 1e-9)
    return Alo, Ahi


NIVELES = [(1.0, 0), (2.30, 2), (4.61, 2), (5.99, 2)]
ROOTS = {lvl: A_roots(lvl) for lvl, _ in NIVELES}
A1m, A1p = ROOTS[1.0]


def eps_de_A(A, rama):
    s = np.sqrt(np.maximum(A, 0.0))
    return ((s if rama == "+" else -s) - 1.0) * QW / COEF_DV


fig, (aa, ab) = plt.subplots(1, 2, figsize=(6.14, 3.1))

# ================================================================ panel (a)
X_MAX = 3.2
aa.plot(Ap, dchi2p, color=C_GE, lw=1.5)
aa.axvline(A_min, color=GRIS, lw=0.8, ls=":")
for lvl, nd in NIVELES:
    Alo, Ahi = ROOTS[lvl]
    aa.axhline(lvl, color="0.55", lw=0.7, ls="--")
    aa.plot([Alo, Ahi], [lvl, lvl], "o", ms=3.4, color=NEGRO, zorder=5)
    aa.text(X_MAX * 0.985, lvl, rf"$\Delta\chi^2={coma(lvl, nd)}$", ha="right", va="bottom", fontsize=6.6, color="0.4")
aa.text(0.80, 0.45, rf"$A_-={coma(A1m,3)}$", fontsize=6.8, ha="right", va="center", color=NEGRO)
aa.text(1.47, 0.45, rf"$A_+={coma(A1p,3)}$", fontsize=6.8, ha="left", va="center", color=NEGRO)
aa.set_xlim(0, X_MAX); aa.set_ylim(0, 7.5)
aa.set_xlabel(r"$A_{\rm amp}$")
aa.set_ylabel(r"$\Delta\chi^2(A_{\rm amp})$ (perfilando $\alpha$)")
aa.set_title("(a) perfil en amplitud", fontsize=9.5)
formatear_coma(aa)

# ================================================================ panel (b)
def fmt_pm(x, nd=3):
    s = f"{abs(x):.{nd}f}".replace(".", ",")
    return f"{'+' if x >= 0 else '−'}{s}"


eg = np.linspace(-0.12, 0.48, 600)
Ag = (1.0 + COEF_DV * eg / QW) ** 2
ab.plot(eg, Ag, color=C_GE, lw=1.5)
ab.axhspan(A1m, A1p, color=GRIS, alpha=0.18, zorder=0)
ab.text(0.25, 1.1, r"$\Delta\chi^2\leq1$", ha="center", va="center", fontsize=6.8, color="0.35")

ramas = sorted([eps_de_A(A1m, "+"), eps_de_A(A1p, "+")]), sorted([eps_de_A(A1p, "-"), eps_de_A(A1m, "-")])
colores_rama = [C_GE, "#8a4a2f"]
for (e_lo, e_hi), col in zip(ramas, colores_rama):
    ab.axvspan(e_lo, e_hi, color=col, alpha=0.22, zorder=0)
    for e in (e_lo, e_hi):
        ab.plot([e, e], [0.0, (1.0 + COEF_DV * e / QW) ** 2], color=col, lw=0.6, ls=":", zorder=1)
    ab.text(0.5 * (e_lo + e_hi), 5.5, f"[{fmt_pm(e_lo)}, {fmt_pm(e_hi)}]",
            ha="center", va="bottom", fontsize=6.4, color=col)

ab.axhline(1.0, color=GRIS, lw=0.7, ls=":")
ab.plot(0, 1.0, "*", ms=8, color=NEGRO, zorder=6)
nota(ab, (0, 1.0), "SM", (0.07, 0.35), color=NEGRO, fs=7.0)
ab.set_xlim(eg[0], eg[-1]); ab.set_ylim(0, 6.5)
ab.set_xlabel(r"$\varepsilon_{ee}^{dV}$ (Ge)")
ab.set_ylabel(r"$A_{\rm amp}(\varepsilon)$")
ab.set_title("(b) traslado al plano NSI", fontsize=9.5)
formatear_coma(ab)

fig.tight_layout()
out = f"{BASE}/fig_conus_A_eps.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = f"{DEST}/fig_conus_A_eps.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print(f"  A_min={A_min:.4f}  A1-={A1m:.4f}  A1+={A1p:.4f}")
for (e_lo, e_hi), (p_lo, p_hi) in zip(ramas, PUBLICADO):
    print(f"  rama propia=[{e_lo:.4f},{e_hi:.4f}]  publicado=[{p_lo:.4f},{p_hi:.4f}]")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
