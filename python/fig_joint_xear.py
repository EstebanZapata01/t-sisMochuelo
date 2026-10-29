#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combinacion Xe+Ar con exposicion total fija repartida (Cap. 6,
sec:combinacion_xear): E_Ar=f_Ar E_tot, E_Xe=(1-f_Ar) E_tot, E_tot=192 kg dia
(la exposicion real base de Xe). MISMO marco que tab:ideal para los dos
blancos: umbral N_e>=1 (sin corte superior), eps_ROI=1, respuesta central,
SIN fondo (Asimov):
    chi2_t(eps) = (1-A_t(eps))^2 S2_t(E_t),   S2_t(E_t)=E_t*Sum_{Ne>=1} R_k^SM,
    chi2_Xe+Ar(eps) = chi2_Xe(eps) + chi2_Ar(eps).
CONTROL (verificado en stdout, no en la figura): f_Ar=0 y f_Ar=1 con esta
S2 reproducen tab:ideal a 192 kg dia (chequeo con Dchi2=2.706, 1 g.d.l.).

Contornos 90% C.L. (Dchi2=4.605, 2 g.d.l.): solo Xe, solo Ar, f_Ar=1/2. En
trazo fino, el lugar analitico A_t(eps)=1 de cada blanco (dos rectas o una
elipse por el origen y por la solucion de cambio de signo q_ee=-2 Q_W^t).

Panel (a) (eps_ee^dV,eps_ee^uV): A_t=1 son DOS RECTAS PARALELAS (bandas) por
blanco. Panel (b) (eps_ee^dV,eps_emu^dV): A_t=1 es una elipse (anillo) por
blanco -- el mismo target, no cambia de DIRECCION, cambia de CENTRO y RADIO
al pasar de Xe a Ar (Sec. "anillo ciego": eso ya esta en fig_blind_rings).

Entrada: datos/espectro_Ne_ideal_{Xe,Ar}.dat
Salida : datos/fig_joint_xear.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from estilo_tesis import aplicar, C_XE, C_AR, GRIS, NEGRO, formatear_coma
from nsi_core import TARGET, S2_ideal, A_amp_grid, COEF
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8, "lines.linewidth": 1.15,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
E_TOT = 192.0
K = 1
N_GRID = 4001
S2 = {t: S2_ideal(BASE, t, K) * E_TOT for t in ("Xe", "Ar")}

# ---------------------------------------------------------------- control
for t in ("Xe", "Ar"):
    a90 = 1 + np.sqrt(2.706 / S2[t])
    print(f"  control tab:ideal {t} N_e>=1 @192kgd: A90={a90:.4f}")


def chi2_joint2(ex, ey, f_Ar, fx, fy, lx, ly):
    c = np.zeros_like(ex)
    for tag, f in (("Xe", 1 - f_Ar), ("Ar", f_Ar)):
        if f <= 0:
            continue
        Z, N, QW = TARGET[tag]["Z"], TARGET[tag]["N"], TARGET[tag]["QW"]
        A = A_amp_grid(ex, ey, Z, N, QW, fx, fy, lx, ly)
        c = c + (1 - A) ** 2 * (f * S2[tag])
    return c


def lineas_A1(fx, fy, lx, ly, tag, lim):
    Z, N, QW = TARGET[tag]["Z"], TARGET[tag]["N"], TARGET[tag]["QW"]
    cx, cy = COEF[fx](Z, N), COEF[fy](Z, N)
    if lx == ly:
        t = np.array(lim) * 1.5
        return [(t, -cx / cy * t), (t, (-2 * QW - cx * t) / cy)]
    else:
        th = np.linspace(0, 2 * np.pi, 400)
        return [((QW * np.cos(th) - QW) / cx, (QW * np.sin(th)) / cy)]


PANELES = [
    (r"$\varepsilon_{ee}^{dV}$", r"$\varepsilon_{ee}^{uV}$", "dV", "uV", "ee", "ee",
     (-1.0, 1.0), "(a) rectas"),
    (r"$\varepsilon_{ee}^{dV}$", r"$\varepsilon_{e\mu}^{dV}$", "dV", "dV", "ee", "emu",
     (-0.5, 0.5), "(b) circunferencias"),
]
CURVAS = [(0.0, "-", C_XE, "solo Xe"), (1.0, "-", C_AR, "solo Ar"), (0.5, "-", NEGRO, "$f_{\\rm Ar}=1/2$")]

fig, axs = plt.subplots(1, 2, figsize=(5.53, 2.85))
for ax, (xl, yl, fx, fy, lx, ly, lim, tit) in zip(axs, PANELES):
    g = np.linspace(lim[0], lim[1], N_GRID)
    EX, EY = np.meshgrid(g, g)
    for f_Ar, ls, col, lab in CURVAS:
        lw = 2.2 if f_Ar == 0.5 else 1.5
        Z2 = chi2_joint2(EX, EY, f_Ar, fx, fy, lx, ly)
        ax.contour(EX, EY, Z2, levels=[4.605], colors=[col], linestyles=[ls], linewidths=lw)
        ax.plot([], [], color=col, ls=ls, lw=lw, label=lab)
    for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
        for (cxv, cyv) in lineas_A1(fx, fy, lx, ly, tag, lim):
            ax.plot(cxv, cyv, color=col, lw=0.6, alpha=0.55)
    ax.plot(0, 0, "*", ms=13, color=NEGRO, zorder=5)
    ax.axhline(0, color="0.85", lw=0.5); ax.axvline(0, color="0.85", lw=0.5)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(tit, fontsize=9.5)
    formatear_coma(ax)
axs[0].plot([], [], "*", ms=8, color=NEGRO, label="SM")
axs[0].legend(loc="upper right", fontsize=7.5)

# ---- inset cerca del origen (ambos paneles) ----
for ax, (xl, yl, fx, fy, lx, ly, lim, tit) in zip(axs, PANELES):
    zlim = 0.10 * lim[1]
    axins = ax.inset_axes([0.03, 0.03, 0.32, 0.32])
    g2 = np.linspace(-zlim, zlim, 1201)
    EX2, EY2 = np.meshgrid(g2, g2)
    for f_Ar, ls, col, lab in CURVAS:
        lw = 1.8 if f_Ar == 0.5 else 1.2
        Z2 = chi2_joint2(EX2, EY2, f_Ar, fx, fy, lx, ly)
        axins.contour(EX2, EY2, Z2, levels=[4.605], colors=[col], linestyles=[ls], linewidths=lw)
    axins.plot(0, 0, "*", ms=8, color=NEGRO)
    axins.set_xlim(-zlim, zlim); axins.set_ylim(-zlim, zlim); axins.set_aspect("equal")
    axins.set_xticks([]); axins.set_yticks([])
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5", lw=0.6)

# ---- "segunda isla" (panel a): cruce de las 2 rectas de cancelacion, marcado con flecha ----
QWx, QWa = TARGET["Xe"]["QW"], TARGET["Ar"]["QW"]
Amat = np.array([[COEF["dV"](TARGET["Xe"]["Z"], TARGET["Xe"]["N"]), COEF["uV"](TARGET["Xe"]["Z"], TARGET["Xe"]["N"])],
                 [COEF["dV"](TARGET["Ar"]["Z"], TARGET["Ar"]["N"]), COEF["uV"](TARGET["Ar"]["Z"], TARGET["Ar"]["N"])]])
bvec = np.array([-2 * QWx, -2 * QWa])
cx0, cy0 = np.linalg.solve(Amat, bvec)
print(f"  cruce lineas de cancelacion (panel a) = ({cx0:.4f}, {cy0:.4f})")
ax0 = axs[0]
ax0.annotate("isla", xy=(cx0, cy0), xycoords="data",
             xytext=(0.97, -0.02), textcoords="data",
             fontsize=7.5, ha="right", va="bottom", color=NEGRO,
             arrowprops=dict(arrowstyle="-", lw=0.7, color=NEGRO, shrinkB=2,
                              connectionstyle="arc3,rad=-0.15"))

fig.tight_layout()
out = f"{BASE}/fig_joint_xear.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
dest = f"{DEST}/fig_joint_xear.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"  malla: {N_GRID}x{N_GRID} en cada panel; insets 1201x1201")
print(f"  S2 (E_tot=192): Xe={S2['Xe']:.4f}  Ar={S2['Ar']:.4f}")
print(f"  para el pie (se quito la linea de condiciones de la figura):")
print(f"    E_tot={E_TOT:.0f} kg dia (E_Xe=(1-f_Ar)E_tot, E_Ar=f_Ar E_tot);"
      f" umbral N_e>=1, sin corte superior, eff_ROI=1, sin fondo (Asimov)")
print(f"  -> {dest}  (md5 {md5(out)})")
