#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
El umbral en electrones se vuelve un umbral en energia del antineutrino.
Fusiona fig_electrones_medios.py (a) y fig_energia_efectiva.py (b), ancho
completo de pagina.

(a) <N_e^ext> = T Qy(T) EEE de Xe y Ar frente a T (log), horizontales en 1
    y 4 electrones; cruces T1, T4 de cada blanco marcados con su valor.
(b) Distribuciones normalizadas en E_nu (flujo, Ar N_e>=1, Xe N_e>=1, Xe
    N_e>=4), zona >8 MeV sombreada; verticales finas en
    E_min(T_k) = [T_k + sqrt(T_k^2 + 2 M T_k)]/2 para cada cruce de (a),
    del color del blanco (la de Xe k=4 cae fuera del espectro, ~10,3 MeV).

Entrada: datos/nest_{,Ar_}218V_dense.txt, datos/dsigma_bare_{Xe,Ar}.dat,
         datos/energia_efectiva_{Xe,Ar}.dat (fig_energia_efectiva.py)
Salida : datos/fig_umbral_electrones.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, NEGRO, GRIS, coma, formatear_coma
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8, "lines.linewidth": 1.15,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
EEE = {"Xe": 0.328, "Ar": 0.99}
ARCHIVO = {"Xe": "nest_218V_dense.txt", "Ar": "nest_Ar_218V_dense.txt"}
COLOR = {"Xe": C_XE, "Ar": C_AR}
M_MEV = {"Xe": 131.293 * 931.49410242, "Ar": 39.948 * 931.49410242}
T_MIN_PIPELINE = 0.2

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.14, 3.1))

# ================================================================ panel (a)
tmax8 = {}
for tag in ("Xe", "Ar"):
    d = np.loadtxt(f"{BASE}/dsigma_bare_{tag}.dat", comments="#")
    tmax8[tag] = d[d[:, 0] == 8.0, 4][0]

datos_T, datos_Qy = {}, {}
for tag in ("Xe", "Ar"):
    d = np.loadtxt(f"{BASE}/{ARCHIVO[tag]}", comments="#")
    datos_T[tag], datos_Qy[tag] = d[:, 0], d[:, 1]


def cruce(T, Ne, n):
    idx = np.where(Ne >= n)[0]
    if not len(idx):
        return np.nan
    i = idx[0]
    return T[0] if i == 0 else np.interp(n, [Ne[i - 1], Ne[i]], [T[i - 1], T[i]])


T_cruce = {}
for tag in ("Xe", "Ar"):
    T, Qy = datos_T[tag], datos_Qy[tag]
    Ne = T * Qy * EEE[tag]
    m = Ne > 0
    a1.plot(T[m], Ne[m], color=COLOR[tag], lw=1.4, label=tag)
    T_cruce[tag] = {n: cruce(T[m], Ne[m], n) for n in (1, 4)}
    a1.axvline(tmax8[tag], color=COLOR[tag], lw=0.8, ls=":", alpha=0.7)

for n in (1, 4):
    a1.axhline(n, color=NEGRO, lw=0.8, ls="--", alpha=0.6)
a1.text(3.0, 1.3, "$N_e=1$", ha="center", fontsize=7, color=NEGRO)
a1.text(0.24, 4.6, "$N_e=4$", fontsize=7, color=NEGRO)

OFF = {("Xe", 1): (7, -11), ("Xe", 4): (6, 6), ("Ar", 4): (6, -14)}
VA = {("Xe", 1): "top", ("Xe", 4): "bottom", ("Ar", 4): "top"}
HA = {("Xe", 1): "left", ("Xe", 4): "left", ("Ar", 4): "left"}
for tag in ("Xe", "Ar"):
    for n in (1, 4):
        Tk = T_cruce[tag][n]
        a1.plot(Tk, n, "o", ms=5, color=COLOR[tag], zorder=6)
        if (tag, n) == ("Ar", 1):
            a1.text(0.205, 1.25, rf"${coma(Tk,3)}$", ha="left", va="bottom",
                    fontsize=6.6, color=COLOR[tag])
        else:
            a1.annotate(rf"${coma(Tk,3)}$", (Tk, n), textcoords="offset points",
                        xytext=OFF[(tag, n)], ha=HA[(tag, n)], va=VA[(tag, n)],
                        fontsize=6.6, color=COLOR[tag])

a1.text(tmax8["Xe"] * 1.1, 7.0, r"$T^{\mathrm{máx}}(8\,{\rm MeV})$", rotation=90,
        ha="left", va="center", fontsize=6.2, color=C_XE)
a1.text(tmax8["Ar"] * 1.06, 3.0, r"$T^{\mathrm{máx}}(8\,{\rm MeV})$", rotation=90,
        ha="left", va="center", fontsize=6.2, color=C_AR)

a1.set_xscale("log")
a1.set_xlim(T_MIN_PIPELINE, 5.0)
a1.set_ylim(0, 9)
a1.set_xticks([0.2, 0.5, 1, 2, 5]); a1.set_xticks([], minor=True)
a1.set_xlabel(r"Energía de retroceso nuclear, $T$  [keV]")
a1.set_ylabel(r"$\langle N_e^{\rm ext}\rangle = T\,Q_y(T)\,\mathrm{EEE}$")
a1.legend(loc="upper left", fontsize=7.8)
a1.set_title("(a) umbral en electrones", fontsize=9.5)
formatear_coma(a1)

# ================================================================ panel (b)
d_xe = np.loadtxt(f"{BASE}/energia_efectiva_Xe.dat", comments="#")
d_ar = np.loadtxt(f"{BASE}/energia_efectiva_Ar.dat", comments="#")
E_xe, flujo, xe_ne1, xe_ne4 = d_xe[:, 0], d_xe[:, 1], d_xe[:, 2], d_xe[:, 3]
E_ar, ar_ne1, ar_ne4 = d_ar[:, 0], d_ar[:, 2], d_ar[:, 3]


def normalizar(E, y):
    return y / np.trapz(y, E)


CURVAS = [
    ("Flujo del reactor", E_xe, normalizar(E_xe, flujo), NEGRO, "-"),
    (r"Xe, $N_e\geq1$", E_xe, normalizar(E_xe, xe_ne1), C_XE, "-"),
    (r"Xe, $N_e\geq4$", E_xe, normalizar(E_xe, xe_ne4), C_XE, "--"),
    (r"Ar, $N_e\geq1$", E_ar, normalizar(E_ar, ar_ne1), C_AR, "-"),
    (r"Ar, $N_e\geq4$", E_ar, normalizar(E_ar, ar_ne4), C_AR, "--"),
]

a2.axvspan(8.0, 10.0, color=GRIS, alpha=0.15, zorder=0)
a2.text(9.0, 0.55, "forma\nextrap.", transform=a2.get_xaxis_transform(),
        ha="center", va="top", fontsize=6.6, color="0.35")
for lab, E, y, col, ls in CURVAS:
    a2.plot(E, y, color=col, ls=ls, lw=1.2, label=lab)
a2.set_xlim(0, 10); a2.set_ylim(bottom=0)

Emin = {}
for tag in ("Xe", "Ar"):
    for n in (1, 4):
        Tk_MeV = T_cruce[tag][n] / 1000.0
        Emin[(tag, n)] = (Tk_MeV + np.sqrt(Tk_MeV**2 + 2 * M_MEV[tag] * Tk_MeV)) / 2
        e = Emin[(tag, n)]
        if e <= 10.0:
            a2.plot([e, e], [0.0, 0.42], color=COLOR[tag], lw=0.8, ls=(0, (4, 3)), alpha=0.8)
a2.annotate(rf"$E_\nu^{{\mathrm{{mín}}}}(T_4^{{\rm Xe}})={coma(Emin[('Xe',4)],1)}$ MeV",
            xy=(10.0, 0.93), xycoords=("data", "axes fraction"),
            xytext=(9.3, 0.93), textcoords=("data", "axes fraction"),
            ha="right", va="center", fontsize=6.2, color=C_XE,
            arrowprops=dict(arrowstyle="->", lw=0.7, color=C_XE, shrinkA=0, shrinkB=2))

a2.set_xlabel(r"Energía del antineutrino, $E_\nu$  [MeV]")
a2.set_ylabel(r"Distribución normalizada  [MeV$^{-1}$]")
a2.legend(loc="center", bbox_to_anchor=(0.33, 0.68), fontsize=6.6, framealpha=0.9)
a2.set_title("(b) umbral en $E_\\nu$", fontsize=9.5)
formatear_coma(a2)

fig.tight_layout()
out = f"{BASE}/fig_umbral_electrones.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = f"{DEST}/fig_umbral_electrones.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print("  T cruce [keV] y E_min [MeV]:")
for tag in ("Xe", "Ar"):
    for n in (1, 4):
        print(f"    {tag} k={n}: T={T_cruce[tag][n]:.4f} keV  E_min={Emin[(tag,n)]:.4f} MeV")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
