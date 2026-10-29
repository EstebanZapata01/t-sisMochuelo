#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
El fondo fija sigma_i, y la senal del SM queda 2-3 ordenes por debajo.
Reemplaza fig_xe_senal_fondo.py y fig8_residuo_ON_OFF.pdf en dos paneles
con el mismo eje x (energia corregida, 110-190 PE).

Arriba: residuo ON-OFF observado (chi2.f90) con la banda de exclusion al
90% C.L. (igual que chi2RED100.py).
Abajo, log-y: fondo B_i reconstruido de la muestra OFF (misma reconstruccion
de sensibilidad_exposicion_completa.py: R_ON-R_OFF=delta_ON_OFF,
R_ON/T_ON0+R_OFF/T_OFF0=sigma_stat^2), sigma_i del residuo y la senal R_i
del SM (chi2.f90).

Entrada: datos/chi2_ON_OFF_banda.dat (PE_center, R_SM, lim90, A90*R_SM,
         delta_ON_OFF, sigma_stat; cabecera A_best, A_90)
Salida : datos/fig_xe_residuo_fondo.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
from estilo_tesis import aplicar, C_XE, NEGRO, GRIS, formatear_coma
from leer_fortran import leer_dat
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8, "lines.linewidth": 1.15,
})

DATADIR = Path("/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos")
DEST = Path("/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia")
T_ON0, T_OFF0 = 192.0, 61.0   # kg dia, exposiciones FIDUCIARIAS 2024

banda, meta = leer_dat(DATADIR / "chi2_ON_OFF_banda.dat")
pe, R_SM, lim, dN, sig = banda[:, 0], banda[:, 1], banda[:, 2], banda[:, 4], banda[:, 5]
ancho = (189.0 - 110.0) / len(pe)

# reconstruccion del fondo OFF (misma cuenta que sensibilidad_exposicion_completa.py)
R_OFF = (sig**2 - dN / T_ON0) / (1.0 / T_ON0 + 1.0 / T_OFF0)
n_neg = int((R_OFF < 0).sum())
R_OFF = np.maximum(R_OFF, 0.0)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(5.22, 4.55), sharex=True,
                              gridspec_kw=dict(height_ratios=[1.3, 1.0], hspace=0.06))

a1.bar(pe, lim, width=ancho, color=C_XE, alpha=0.35, lw=0, label=r"límite al $90\%$ C.L.")
a1.errorbar(pe, dN, yerr=sig, fmt="o", color=NEGRO, ms=4, lw=0.9, capsize=2.0, label="ON$-$OFF")
a1.axhline(0, color=GRIS, lw=0.7)
a1.set_ylim(-0.5, 0.5)
a1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
a1.set_ylabel(r"eventos $\cdot$ kg$^{-1}\cdot$ día$^{-1}$")
a1.legend(loc="upper right", fontsize=7.8)

assert (R_OFF > sig).all()   # la barra +-sigma_i cabe en la escala log
bordes = np.linspace(110.0, 189.0, 16)
assert np.allclose(0.5 * (bordes[1:] + bordes[:-1]), pe, atol=0.01)
a2.stairs(R_OFF, bordes, baseline=None, color="0.35", lw=1.1,
          label=r"fondo $B_i$ (reconstruido del OFF), barra $\pm\sigma_i$ del residuo")
a2.errorbar(pe, R_OFF, yerr=sig, fmt="none", ecolor=GRIS, elinewidth=0.9, capsize=1.8, capthick=0.9)
a2.stairs(R_SM, bordes, baseline=None, color=C_XE, lw=1.3, label=r"$R_{\rm SM}$ (predicción propia)")
a2.set_yscale("log")
a2.set_xlim(108, 191)
a2.xaxis.set_major_locator(ticker.MultipleLocator(10))
a2.set_xlabel("energía corregida [PE]")
a2.set_ylabel(r"eventos $\cdot$ kg$^{-1}\cdot$ día$^{-1}$")
a2.legend(loc="center left", bbox_to_anchor=(0.01, 0.49), fontsize=7.5)
formatear_coma(a1, x=False); formatear_coma(a2, x=False)

fig.tight_layout()
out = DATADIR / "fig_xe_residuo_fondo.png"
fig.savefig(out, dpi=300)
plt.close(fig)
dest = DEST / "fig_xe_residuo_fondo.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print(f"  bines con R_OFF reconstruido negativo (capados a 0): {n_neg}/{len(pe)}")
print(f"  B_i rango: {R_OFF.min():.4f} - {R_OFF.max():.4f}   R_SM rango: {R_SM.min():.2e} - {R_SM.max():.2e}")
print(f"  cociente tipico B_i/R_SM: ~{np.median(R_OFF[R_SM>0]/R_SM[R_SM>0]):.0f}x")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
