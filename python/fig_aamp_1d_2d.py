#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lo que en 1D son dos puntos, en 2D es una circunferencia. Reemplaza
fig_aamp_1d.py y fig_misma_prediccion_plano.py en un solo par de paneles.
Dibujada ya al tamano final de pagina (ancho completo, 15,6 cm = 6,14 in)
para que LaTeX no reescale letra.

(a) A_amp(eps_ee^dV) con los demas eps nulos, Xe y Ar: A_amp=1 (SM y la
    segunda solucion 2 eps_0) y A_amp=4, con los dos cortes del Xe
    marcados (eps_ee^dV = -0,179 y +0,538). P1 = el de -0,179.
(b) Plano (eps_ee^dV, eps_emu^dV) del Xe: las mismas curvas de nivel
    A_amp=0,25/1/2/4; P1 (eps_emu=0, el MISMO punto y marcador de (a)),
    el punto gemelo en angulo 180 (mismo marcador sin rotulo que en (a)),
    P2, P3 (A_amp=4, angulos +-55) y P4 (A_amp=1, angulo 180).

Entrada: nsi_core (analitico)
Salida : datos/fig_aamp_1d_2d.png (+ copia en Proyecto_LaTeX/figuras_metodologia/)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, NEGRO, CICLO, coma, formatear_coma, nota
from nsi_core import TARGET, A_amp_grid, punto_en_circulo, COEF
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8, "lines.linewidth": 1.15,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.14, 3.35))

# ================================================================ panel (a)
eps = np.linspace(-0.25, 0.65, 900)
eps0 = {}
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    Z, N, QW = TARGET[tag]["Z"], TARGET[tag]["N"], TARGET[tag]["QW"]
    A = A_amp_grid(eps, 0.0, Z, N, QW, "dV", "dV")
    a1.plot(eps, A, color=col, lw=1.2, label=tag)
    eps0[tag] = -QW / (Z + 2 * N)

a1.axhline(1.0, color=NEGRO, lw=0.8, ls="--")
a1.axhline(4.0, color=NEGRO, lw=0.8, ls=":")
a1.text(0.57, 1.15, "$A_{\\rm amp}=1$", fontsize=7.5, color=NEGRO, ha="right")
a1.text(0.2, 4.15, "$A_{\\rm amp}=4$", fontsize=7.5, color=NEGRO, ha="center", va="bottom")

Z, N, QW = TARGET["Xe"]["Z"], TARGET["Xe"]["N"], TARGET["Xe"]["QW"]
cx = COEF["dV"](Z, N)
x_m, x_p = QW * (2 - 1) / cx, QW * (-2 - 1) / cx   # -0.179 (P1), +0.538 (gemelo)
a1.plot(0, 1, "*", ms=9, color=NEGRO, zorder=6)
a1.plot(eps0["Xe"], 0, "x", ms=7, mew=1.6, color=C_XE, zorder=6)
a1.plot(2 * eps0["Xe"], 1, "D", ms=5.5, color=C_XE, zorder=6)
a1.plot(x_m, 4, "o", ms=6.5, color=C_XE, zorder=6)                                  # P1
a1.plot(x_p, 4, "d", ms=6.5, mfc="white", mec=C_XE, mew=1.4, zorder=6)              # gemelo, sin rotulo
nota(a1, (0, 1), "SM", (0.05, 2.1), color=NEGRO, fs=7.5)
nota(a1, (eps0["Xe"], 0), rf"$\varepsilon_0={coma(eps0['Xe'],4)}$", (-0.22, 0.55), color=C_XE, fs=7.5)
nota(a1, (2 * eps0["Xe"], 1), rf"$2\varepsilon_0={coma(2*eps0['Xe'],3)}$", (0.33, 0.35), color=C_XE, fs=7.5)
nota(a1, (x_m, 4), rf"P$_1$: ${coma(x_m,3)}$", (-0.12, 3.6), color=C_XE, fs=7.5)
a1.set_xlim(-0.25, 0.65); a1.set_ylim(-0.3, 5.6)
a1.set_xlabel(r"$\varepsilon_{ee}^{dV}$"); a1.set_ylabel(r"$A_{\rm amp}$")
a1.legend(loc="upper left", fontsize=7.8)
a1.set_title("(a) un parámetro", fontsize=9.5)
formatear_coma(a1)

# ================================================================ panel (b)
A_STAR = 4.0
Zx, Nx, QWx = TARGET["Xe"]["Z"], TARGET["Xe"]["N"], TARGET["Xe"]["QW"]
eps0x = -QWx / (Zx + 2 * Nx)

PUNTOS = [(r"P$_1$", 0.0, A_STAR), (r"P$_2$", -55.0, A_STAR), (r"P$_3$", 55.0, A_STAR),
          (r"P$_4$", 180.0, 1.0)]
COL4 = [C_XE] + list(CICLO[1:3]) + [NEGRO]
MARK = ["o", "s", "^", "D"]
OFF = [(-0.42, -0.22), (0.20, 0.48), (0.20, -0.48), (0.36, -0.22)]
HAL = ["center", "left", "left", "center"]

g = np.linspace(-0.6, 0.6, 700)
EX, EY = np.meshgrid(g, g)
Amap = A_amp_grid(EX, EY, Zx, Nx, QWx, "dV", "dV")
a2.contour(EX, EY, Amap, levels=[0.25, 2.0], colors=["0.6"], linewidths=0.8)
a2.contour(EX, EY, Amap, levels=[4.0], colors=[CICLO[0]], linewidths=1.2)
a2.contour(EX, EY, Amap, levels=[1.0], colors=[NEGRO], linewidths=1.1)
a2.plot([], [], color="0.6", lw=0.8, label=r"$A_{\rm amp}=0{,}25\ \mathrm{y}\ 2$")
a2.plot([], [], color=NEGRO, lw=1.1, label=r"$A_{\rm amp}=1$")
a2.plot([], [], color=CICLO[0], lw=1.2, label=r"$A_{\rm amp}=4$")
a2.plot(0, 0, "*", ms=9, color=NEGRO, zorder=6)
a2.plot(eps0x, 0, "x", ms=7, mew=1.6, color=NEGRO, zorder=6)
nota(a2, (0, 0), "SM", (-0.16, 0.20), color=NEGRO, fs=7.5, ha="center")

xp_b, yp_b = punto_en_circulo(4.0, QWx, Zx, Nx, 180.0)
a2.plot(xp_b, yp_b, "d", ms=6.5, mfc="white", mec=C_XE, mew=1.4, zorder=6)   # gemelo de P1, sin rotulo

for (lab, ang, A_t), col, mk, off, ha in zip(PUNTOS, COL4, MARK, OFF, HAL):
    px, py = punto_en_circulo(A_t, QWx, Zx, Nx, ang)
    a2.plot(px, py, mk, ms=6.5, color=col, zorder=6)
    nota(a2, (px, py), lab, off, color=col, fs=7.8, ha=ha)
a2.set_xlim(-0.6, 0.6); a2.set_ylim(-0.6, 0.6); a2.set_aspect("equal")
a2.set_xlabel(r"$\varepsilon_{ee}^{dV}$"); a2.set_ylabel(r"$\varepsilon_{e\mu}^{dV}$")
a2.legend(loc="upper left", fontsize=7.5)
a2.set_title("(b) dos parámetros (Xe)", fontsize=9.5)
formatear_coma(a2)

fig.tight_layout()
out = f"{BASE}/fig_aamp_1d_2d.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = f"{DEST}/fig_aamp_1d_2d.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print(f"  eps_0(Xe)={eps0['Xe']:.4f}  eps_0(Ar)={eps0['Ar']:.4f}")
print(f"  Xe, cortes con A_amp=4: {x_m:.4f} (P1, angulo 0) y {x_p:.4f} (angulo 180)")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
