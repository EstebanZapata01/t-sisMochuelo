#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kappa_eee_plot.py -- figura del factor de amplificacion kappa_EEE.

SOLO GRAFICACION. El calculo vive en python/kappa_eee.py (hay que correrlo
antes: genera datos/eee_kappa_grid_{Xe,Ar}.dat y datos/kappa_eee.dat).

Un solo panel: kappa_EEE vs. umbral del ROI n_thr, para Xe y Ar. Linea =
metodo analitico (identidad binomial exacta); anillos abiertos = diferencia
finita sobre el pipeline compilado (coinciden a <0.02%). Estrellas = piso
fisico de cada blanco (Xe n_thr=4, Ar n_thr=1).

Salida: datos/fig_kappa_eee.png
"""
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar
aplicar()
import matplotlib.ticker as ticker

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kappa_eee import kappa_analytic, EEE_NOM, EEE_SIGREL, NTHR_PHYS

DATOS = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
C = {"Xe": "#33546e", "Ar": "#a86a43"}
NTHRS = (1, 2, 3, 4)


def read_fd_kappa():
    """{tag: {n_thr: kappa_dif_finita}} desde datos/kappa_eee.dat."""
    txt = open(f"{DATOS}/kappa_eee.dat", encoding="utf-8").read()
    out = {}
    for tag in ("Xe", "Ar"):
        m = re.search(rf"\[{tag}\] ideal.*?(?=\n\n|\n  \[)", txt, re.S)
        blk = m.group(0) if m else ""
        out[tag] = {int(a): float(b) for a, b in
                    re.findall(r"n_thr=(\d):\s*kappa_analitico\s*=\s*[-0-9.]+\s*"
                               r"kappa_dif\.finita\s*=\s*([-0-9.]+)", blk)}
    return out


def main():
    fd = read_fd_kappa()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))

    ax.axhline(1.0, color="0.6", lw=1.0, ls="--")
    ax.text(1.05, 1.03, r"normalización pura ($\kappa=1$)", fontsize=8,
            color="0.45", va="bottom")

    for tag in ("Xe", "Ar"):
        ka = [kappa_analytic(tag, EEE_NOM[tag], n)[0] for n in NTHRS]
        ax.plot(NTHRS, ka, "-o", color=C[tag], ms=5, label=f"{tag} (analítico)")
        kf = [fd[tag].get(n, np.nan) for n in NTHRS]
        ax.plot(NTHRS, kf, "o", mfc="none", mec=C[tag], ms=10, mew=1.2,
                label=f"{tag} (dif. finita)")
        nphys = NTHR_PHYS[tag]
        ax.plot(nphys, kappa_analytic(tag, EEE_NOM[tag], nphys)[0], marker="*",
                color=C[tag], ms=17, mec="black", mew=0.6, zorder=5)

    kx4 = kappa_analytic("Xe", EEE_NOM["Xe"], 4)[0]
    ka1 = kappa_analytic("Ar", EEE_NOM["Ar"], 1)[0]
    ax.annotate(rf"$n_{{\rm thr}}=4$ (Xe): $\kappa={kx4:.2f}$" "\n"
                rf"$\kappa\,\sigma_{{\rm EEE}}/\mathrm{{EEE}}\approx{kx4*EEE_SIGREL['Xe']*100:.0f}\%$",
                xy=(4, kx4), xytext=(1.7, 3.0), fontsize=8.6, color=C["Xe"],
                arrowprops=dict(arrowstyle="->", color=C["Xe"], lw=0.8))
    ax.annotate(rf"$n_{{\rm thr}}=1$ (Ar): $\kappa={ka1:.2f}$",
                xy=(1, ka1), xytext=(1.9, 0.30), fontsize=8.6, color=C["Ar"],
                ha="center",
                arrowprops=dict(arrowstyle="->", color=C["Ar"], lw=0.8))

    ax.set_xlabel(r"Umbral del ROI  $n_{\rm thr}$  ($N_e \geq n_{\rm thr}$)")
    ax.set_ylabel(r"$\kappa_{\rm EEE} = (\mathrm{EEE}/S)\;dS/d\mathrm{EEE}$")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.set_xlim(0.8, 4.4)
    ax.set_ylim(0, 3.8)
    ax.legend(loc="upper left", ncol=2, columnspacing=1.0, handletextpad=0.4)

    fig.tight_layout()
    out = f"{DATOS}/fig_kappa_eee.png"
    fig.savefig(out)
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
