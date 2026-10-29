#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rendimiento de carga Q_y(T) y factor de Fano F(T) de Xe (NEST) y Ar
(LArNEST reescalado a ReD), en dos paneles apilados con el mismo eje T
(log, 0,1-10 keV). Fusiona fig_nest_qy.png y fig_nest_fano.png, 0,8 del
ancho de pagina.

Ar: el tramo con dato medido de ReD 2025 (energias medias, 2,4-7,6 keV) va
en continuo; por debajo y por encima es extrapolacion de LArNEST, discontinuo.

Barras: rango de retrocesos del CEvNS de cada blanco, de 0,2 keV (piso de
la tabla) a su T_max propio a E_nu=8 MeV (borde duro del espectro del
reactor; dsigma_bare_XeAr.py).

Xe: por debajo de T=0,2088 keV la tabla de NEST da Qy=0 (sin ionizacion)
y por convencion F=1 como valor por defecto (relacion Var(Ne)/<Ne>
indefinida con <Ne>=0); se excluye del panel de Fano con el mismo filtro
Qy>0 que ya usa el panel de Qy, para no mostrar ese salto espurio 1->0,5.

Entradas: datos/nest_218V_dense.txt, datos/nest_Ar_218V_dense.txt,
          datos/dsigma_bare_{Xe,Ar}.dat (T_max a 8 MeV)
Salida : datos/fig_nest_qy_fano.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, GRIS, formatear_coma
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8, "lines.linewidth": 1.15,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
T_MIN, T_MAX = 0.1, 10.0
T_MED_MIN, T_MED_MAX = 2.4, 7.6   # keV, rango medido por ReD 2025 (Ar)

dXe = np.loadtxt(f"{BASE}/nest_218V_dense.txt", comments="#")
dAr = np.loadtxt(f"{BASE}/nest_Ar_218V_dense.txt", comments="#")

tmax8 = {}
for tag in ("Xe", "Ar"):
    d = np.loadtxt(f"{BASE}/dsigma_bare_{tag}.dat", comments="#")
    tmax8[tag] = d[d[:, 0] == 8.0, 4][0]

fig, (aq, af) = plt.subplots(2, 1, figsize=(4.91, 4.6), sharex=True,
                              gridspec_kw=dict(height_ratios=[1.2, 1.0], hspace=0.08))

Y_BAR_XE, Y_BAR_AR = 10.0, 9.3
for tag, col, y in (("Xe", C_XE, Y_BAR_XE), ("Ar", C_AR, Y_BAR_AR)):
    aq.annotate("", xy=(0.2, y), xytext=(tmax8[tag], y),
                arrowprops=dict(arrowstyle="|-|", lw=1.1, color=col, mutation_scale=2.8))
    tmax_txt = f"{tmax8[tag]:.2f}".replace(".", ",")
    aq.text(tmax8[tag] * 1.1, y, f"{tag}: 0,2–{tmax_txt} keV",
            ha="left", va="center", fontsize=6.6, color=col)

T_PISO_NEST = 0.2
for a in (aq, af):
    a.axvspan(T_MIN, T_PISO_NEST, color="0.85", alpha=0.6, zorder=0)
    a.axvline(T_PISO_NEST, color="0.4", lw=0.9, ls=":")
aq.text(np.sqrt(T_MIN * T_PISO_NEST), 7.0, "por debajo\ndel mínimo\nde NEST",
        ha="center", va="center", fontsize=6.6, color="0.35")

selXe = (dXe[:, 0] <= T_MAX) & (dXe[:, 0] >= T_MIN) & (dXe[:, 1] > 0)
aq.plot(dXe[selXe, 0], dXe[selXe, 1], color=C_XE, lw=1.4, label="Xe (NEST)")

T_ar, Qy_ar = dAr[:, 0], dAr[:, 1]
sel_panel = (T_ar <= T_MAX) & (T_ar >= T_MIN)
ext = (T_ar < T_MED_MIN) | (T_ar > T_MED_MAX)
aq.plot(T_ar[sel_panel & ext], Qy_ar[sel_panel & ext], color=C_AR, lw=1.2, ls="--",
        label="Ar, extrapolación")
aq.plot(T_ar[sel_panel & ~ext], Qy_ar[sel_panel & ~ext], color=C_AR, lw=1.4,
        label="Ar, rango medido (ReD)")

aq.set_xscale("log")
aq.set_xlim(T_MIN, T_MAX)
aq.set_ylim(0, 10.8)
aq.set_ylabel(r"$Q_y(T)$  [e$^-$/keV]")
aq.legend(loc="lower right", fontsize=7.2)

for d, sel, c, tag in ((dXe, selXe, C_XE, "Xe"), (dAr, sel_panel, C_AR, "Ar")):
    af.plot(d[sel, 0], d[sel, 2], color=c, lw=1.4, label=tag)
af.set_xscale("log")
af.set_xlim(T_MIN, T_MAX); af.set_ylim(0, 1.05)
af.set_xlabel(r"$T$  [keV]")
af.set_ylabel(r"$\mathcal{F}(T)=\mathrm{Var}(N_e)/\langle N_e\rangle$")
af.legend(loc="lower right", fontsize=7.8)

formatear_coma(aq); formatear_coma(af)
fig.tight_layout()
out = f"{BASE}/fig_nest_qy_fano.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = f"{DEST}/fig_nest_qy_fano.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)
print(f"  T_max(8 MeV): Xe={tmax8['Xe']:.3f} keV  Ar={tmax8['Ar']:.3f} keV")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
