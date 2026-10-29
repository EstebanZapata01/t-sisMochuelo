#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Seccion eficaz TOTAL sigma(E_nu) = integral de dsigma/dT entre T=0 y T_max(E_nu)
(FORTRAN90/N_EventosCEvNS_NSI{Xe,Ar}/sigma_total.f90, rutina dsigma_dT real),
frente a la aproximacion lineal sigma_lin=(G_F^2/pi) Q_W^2 E_nu^2 que se cita
en el texto pero no se habia graficado.

Panel izquierdo: sigma(E_nu) por nucleo, Xe y Ar. Helm (linea continua,
el default de produccion), F^2=1 (discontinua), aproximacion lineal
(punteada, point-like, F^2=1 por construccion).
Panel derecho: sigma_Xe/sigma_Ar con Helm frente a E_nu, con la referencia
point-like (Q_W,Xe/Q_W,Ar)^2.

Entrada: datos/sigma_total_{Xe,Ar}.dat (Helm), sigma_total_bare_{Xe,Ar}.dat (F^2=1)
Salida : datos/fig_sigma_total.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, GRIS, formatear_coma, coma
from nsi_core import TARGET
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0, "lines.linewidth": 1.3,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
COLT = {"Xe": C_XE, "Ar": C_AR}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.14, 3.1))

sig = {}
for tag in ("Xe", "Ar"):
    dhelm = np.loadtxt(f"{BASE}/sigma_total_{tag}.dat", comments="#")
    dbare = np.loadtxt(f"{BASE}/sigma_total_bare_{tag}.dat", comments="#")
    sig[tag] = dhelm
    a1.plot(dhelm[:, 0], dhelm[:, 1], color=COLT[tag], lw=1.5, label=f"{tag}, Helm")
    a1.plot(dbare[:, 0], dbare[:, 1], color=COLT[tag], lw=1.2, ls="--", label=rf"{tag}, $F^2=1$")
    a1.plot(dhelm[:, 0], dhelm[:, 2], color=COLT[tag], lw=1.0, ls=":", label=rf"{tag}, aprox. lineal")
a1.set_yscale("log")
a1.set_xlabel(r"$E_\nu$  [MeV]"); a1.set_ylabel(r"$\sigma(E_\nu)$  [cm$^2$]")
a1.set_xlim(0, 10)
a1.legend(loc="lower right", fontsize=7.4)
formatear_coma(a1, y=False)
a1.yaxis.set_major_formatter(mpl.ticker.LogFormatterSciNotation())
a1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ratio = sig["Xe"][:, 1] / sig["Ar"][:, 1]
ref = (TARGET["Xe"]["QW"] / TARGET["Ar"]["QW"]) ** 2
a2.plot(sig["Xe"][:, 0], ratio, color="0.15", lw=1.5, label=r"$\sigma_{\rm Xe}/\sigma_{\rm Ar}$ (Helm)")
a2.axhline(ref, color=GRIS, lw=1.1, ls="--", label=rf"$(Q_{{W,{{\rm Xe}}}}/Q_{{W,{{\rm Ar}}}})^2={coma(ref,1)}$")
a2.set_xlabel(r"$E_\nu$  [MeV]"); a2.set_ylabel(r"$\sigma_{\rm Xe}/\sigma_{\rm Ar}$")
a2.set_xlim(0, 10); a2.set_ylim(12.0, 13.0)
a2.legend(loc="lower right", fontsize=8.0)
formatear_coma(a2)

fig.tight_layout()
out = f"{BASE}/fig_sigma_total.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
dest = f"{DEST}/fig_sigma_total.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

i5 = np.argmin(np.abs(sig["Xe"][:, 0] - 5.0))
j5 = np.argmin(np.abs(sig["Ar"][:, 0] - 5.0))
print(f"  sigma_Xe(5 MeV, Helm) = {sig['Xe'][i5,1]:.4e} cm^2")
print(f"  sigma_Ar(5 MeV, Helm) = {sig['Ar'][j5,1]:.4e} cm^2")
bare = {t: np.loadtxt(f"{BASE}/sigma_total_bare_{t}.dat", comments="#") for t in ("Xe", "Ar")}
for tag in ("Xe", "Ar"):
    dev = np.abs(bare[tag][:, 1] / bare[tag][:, 2] - 1.0)
    print(f"  {tag}: max desviacion |F^2=1/lineal - 1| hasta 10 MeV = {dev.max():.2e} (en E_nu={bare[tag][dev.argmax(),0]:.2f} MeV)")
    red = 1.0 - sig[tag][:, 1] / bare[tag][:, 1]
    print(f"  {tag}: reduccion de sigma_total por Helm (1 - Helm/F^2=1), max hasta 10 MeV = {red.max():.4f} "
          f"(en E_nu={sig[tag][red.argmax(),0]:.2f} MeV); en 10 MeV = {red[-1]:.4f}")
print(f"  ratio Xe/Ar (Helm) en 5 MeV: {ratio[i5]:.3f}  (referencia point-like (QW,Xe/QW,Ar)^2={ref:.3f})")
print(f"  -> {dest}  (md5 {md5(out)})")
