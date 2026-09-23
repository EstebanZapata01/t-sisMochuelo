#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Residuo ON-OFF de RED-100 con la banda de exclusion al 90% C.L. y perfil fisico
(A >= 0) Delta chi^2(A) de la amplitud CEvNS, observado y Asimov con la misma chi2(A), a partir de la salida de chi2.f90.

Entradas: datos/chi2_ON_OFF_banda.dat (PE R_SM A90*R_SM -A90*R_SM dN sigma; cabecera A_best, A_90)
          datos/chi2_ON_OFF_perfil.dat (A chi2(A))
Salidas : datos/fig8_residuo_ON_OFF.pdf, datos/fig9_chi2_perfil.pdf
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from estilo_tesis import aplicar, C_XE, NEGRO, GRIS, FIG15
from leer_fortran import leer_dat
aplicar(grande=True)

DATADIR = Path("/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos")


banda, meta = leer_dat(DATADIR / "chi2_ON_OFF_banda.dat")
perfil, _ = leer_dat(DATADIR / "chi2_ON_OFF_perfil.dat")
pe, R_SM, lim, dN, sig = banda[:, 0], banda[:, 1], banda[:, 2], banda[:, 4], banda[:, 5]
A, chi2 = perfil[:, 0], perfil[:, 1]
A90 = meta["A_90"]
ancho = (189.0 - 110.0) / len(pe)                       # ancho de bin del histograma [PE]

# ---- residuo ON-OFF (energia corregida)
fig, ax = plt.subplots(figsize=FIG15)
ax.bar(pe, lim, width=ancho, color=C_XE, alpha=0.35, lw=0, label=r"límite al $90\%$ C.L.")
ax.errorbar(pe, dN, yerr=sig, fmt="o", color=NEGRO, ms=5.5, lw=1.1, capsize=2.5, label="ON$-$OFF")
ax.axhline(0, color=GRIS, lw=0.7)
ax.set_xlim(108, 191); ax.set_ylim(-0.5, 0.5)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.set_xlabel("energía corregida [PE]")
ax.set_ylabel(r"eventos $\cdot$ kg$^{-1}\cdot$ día$^{-1}$")
ax.legend(loc="upper right")
fig.tight_layout(); fig.savefig(DATADIR / "fig8_residuo_ON_OFF.pdf"); plt.close(fig)

# ---- perfil fisico de Delta chi^2 en la amplitud (A >= 0)
# Observado y Asimov pasan por la MISMA chi2(A) = sum((D - A R_SM)^2/sigma^2), con los mismos R_SM y sigma;
# solo cambian los datos: D = ON-OFF (observado) o D = R_SM (Asimov, A_best = 1).
assert A[0] == 0.0, "el perfil debe partir de A = 0 (A >= 0)"
CL90 = 2.706


def chi2_A(a, D):
    return (((D[None, :] - a[:, None] * R_SM[None, :]) / sig[None, :]) ** 2).sum(axis=1)


a = np.linspace(0, 100, 2001)
d_obs = chi2_A(a, dN) - chi2_A(np.zeros(1), dN)[0]                 # frontera fisica: chi2(A) - chi2(0)
d_asi = chi2_A(a, R_SM) - chi2_A(np.array([1.0]), R_SM)[0]         # minimo en A = 1


def cruce(y):
    """A donde y = 2.706 (interpolacion lineal, tramo creciente)."""
    k = np.argmax(y >= CL90)
    return a[k - 1] + (CL90 - y[k - 1]) * (a[k] - a[k - 1]) / (y[k] - y[k - 1])


A90_obs, A90_esp = cruce(d_obs), cruce(d_asi)
assert abs(A90_obs - A90) < 0.05, (A90_obs, A90)                    # coincide con el A_90 que escribe chi2.f90

fig, ax = plt.subplots(figsize=FIG15)
ax.axhline(CL90, color=NEGRO, lw=1.1, label=r"corte al $90\%$ C.L. ($\Delta\chi^{2}=2{,}71$)")
ax.axhline(1.0, color=NEGRO, lw=0.9, ls=":")
ax.plot(a, d_asi, color=GRIS, ls="-.", label=rf"esperado (Asimov): $A_{{90}}={A90_esp:.0f}$")
ax.plot(a, d_obs, color=C_XE, label=rf"observado (ON$-$OFF): $A_{{90}}={A90_obs:.0f}$")
for x, col in ((A90_esp, GRIS), (A90_obs, C_XE)):
    ax.vlines(x, 0, CL90, color=col, ls=":", lw=1.4)
    ax.plot([x], [CL90], "o", mfc="white", mec=col, mew=1.6, ms=8, zorder=5)
ax.set_xlim(0, 100); ax.set_ylim(0, 5)
ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.set_xlabel(r"amplitud $A$  [$\times$SM]")
ax.set_ylabel(r"$\Delta\chi^{2}$")
ax.legend(loc="upper left")
fig.tight_layout(); fig.savefig(DATADIR / "fig9_chi2_perfil.pdf"); plt.close(fig)

print(f"  A_best fisico={meta.get('A_best', float('nan')):.3f};  A_90 obs={A90_obs:.2f}  A_90 esp={A90_esp:.2f}  "
      f"chi2(0)={chi2_A(np.zeros(1), dN)[0]:.2f}/{len(pe)}")
