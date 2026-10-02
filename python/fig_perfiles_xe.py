#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perfiles Delta chi^2(A) del xenon de RED-100 en una sola figura:
  (1) el motor estadistico alimentado con los tres histogramas y la senal que publica
      la colaboracion (observado y esperado), frente a sus curvas publicadas;
  (2) la cadena propia: histograma de energia y prediccion propia (observado y Asimov).
Funde fig_validacion_red100.png (valida_tabla1.py) y fig9_chi2_perfil.pdf (chi2RED100.py),
con las mismas entradas y las mismas formulas que esos dos scripts.

Entrada: datos/red100_residuo_ONOFF_fig8.csv, datos/red100_perfil_chi2_fig9.csv,
         datos/chi2_ON_OFF_banda.dat (chi2.f90)
Salida : datos/fig_perfiles_xe.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from estilo_tesis import aplicar, C_XE, GRIS, NEGRO, FIG15, formatear_coma
from leer_fortran import leer_dat
aplicar(grande=True)

BASE = "../datos"
CL90 = 2.706
OUT = f"{BASE}/fig_perfiles_xe.png"

# ---- (1) insumos de la colaboracion: tres histogramas, senal publicada R_i = barra_i/63
f = [l.strip().split(",") for l in open(f"{BASE}/red100_residuo_ONOFF_fig8.csv")
     if not l.startswith(("#", "panel"))]
v = np.array([[float(x) for x in r[1:]] for r in f])
D, sig, R = v[:, 1], v[:, 2], v[:, 3] / 63.0
S1, S2, S3 = (D * R / sig**2).sum(), (R**2 / sig**2).sum(), (D**2 / sig**2).sum()
A = np.linspace(0, 130, 5201)
c_obs = S3 - 2 * A * S1 + A**2 * S2
col_obs = c_obs - c_obs.min()                  # minimo fisico en A >= 0
col_esp = S2 * (A - 1) ** 2

pub = {}
for l in open(f"{BASE}/red100_perfil_chi2_fig9.csv"):
    if not l.startswith(("#", "curva")):
        c, a_, d_ = l.strip().split(",")
        pub.setdefault(c, []).append((float(a_), float(d_)))

# ---- (2) cadena propia: histograma de energia, prediccion propia
banda, meta = leer_dat(f"{BASE}/chi2_ON_OFF_banda.dat")
R_SM, dN, sg = banda[:, 1], banda[:, 4], banda[:, 5]
A_BEST = meta["A_best"]


def chi2_A(a, Dd):
    return (((Dd[None, :] - a[:, None] * R_SM[None, :]) / sg[None, :]) ** 2).sum(axis=1)


pro_obs = chi2_A(A, dN) - chi2_A(np.array([A_BEST]), dN)[0]
pro_asi = chi2_A(A, R_SM) - chi2_A(np.array([1.0]), R_SM)[0]


def cruce(y, a_min):
    k = np.argmax((y >= CL90) & (A > a_min))
    return A[k - 1] + (CL90 - y[k - 1]) * (A[k] - A[k - 1]) / (y[k] - y[k - 1])


x_col_esp, x_col_obs = cruce(col_esp, 1.0), cruce(col_obs, S1 / S2)
x_pro_asi, x_pro_obs = cruce(pro_asi, 1.0), cruce(pro_obs, A_BEST)
print(f"3 hist + senal RED-100 : esp {x_col_esp:.1f} obs {x_col_obs:.1f}")
print(f"energia + prediccion   : esp {x_pro_asi:.1f} obs {x_pro_obs:.1f}")

C1 = "#4d4d4d"
fig, ax = plt.subplots(figsize=(FIG15[0], FIG15[1] * 1.08))
ax.axhline(CL90, color=NEGRO, lw=1.0)
ax.text(2.0, CL90 + 0.08, r"$90\%$ C.L.", fontsize=10, va="bottom")
curvas = [(col_esp, col_obs, C1, 1.8, "tres histogramas, señal de RED-100"),
          (pro_asi, pro_obs, C_XE, 2.2, "histograma de energía, predicción propia")]
hs = []
for esp, obs, col, lw, lab in curvas:
    ax.plot(A, esp, color=col, ls="-.", lw=lw)
    h, = ax.plot(A, obs, color=col, ls="-", lw=lw, label=lab)
    hs.append(h)
for k in ("esperada", "observada"):
    q = np.array(pub[k])
    hp, = ax.plot(q[:, 0], q[:, 1], color=NEGRO, ls=":", lw=1.7)
hs.insert(1, hp)
labs = [curvas[0][4], "curvas publicadas por la colaboración", curvas[1][4]]
for x, col in ((x_col_esp, C1), (x_col_obs, C1), (x_pro_asi, C_XE), (x_pro_obs, C_XE)):
    ax.vlines(x, 0, CL90, color=col, ls=":", lw=1.1)
    ax.plot([x], [CL90], "o", mfc="white", mec=col, mew=1.6, ms=7, zorder=5)
h1, = ax.plot([], [], color="0.35", ls="-", lw=1.6)
h2, = ax.plot([], [], color="0.35", ls="-.", lw=1.6)
leg1 = ax.legend(hs, labs, loc="upper left", fontsize=9.4)
ax.add_artist(leg1)
ax.legend([h1, h2], ["observado", "esperado (Asimov)"], loc="lower right", fontsize=9.4)
ax.set_xlim(0, 100); ax.set_ylim(0, 5)
ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.set_xlabel(r"amplitud $A$  [$\times$SM]")
ax.set_ylabel(r"$\Delta\chi^{2}$")
formatear_coma(ax)
fig.tight_layout(); fig.savefig(OUT, dpi=300); plt.close(fig)
print("->", OUT)
