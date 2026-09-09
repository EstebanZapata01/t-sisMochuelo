#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kappa_eee_plot.py -- figura del factor de amplificacion kappa_EEE.

SOLO GRAFICACION. El calculo vive en python/kappa_eee.py (hay que correrlo
antes: genera datos/eee_kappa_grid_{Xe,Ar}.dat y datos/kappa_eee.dat).

Panel A : kappa_EEE vs. umbral del ROI n_thr, para Xe y Ar (metodo analitico
          como linea; diferencia finita como marcadores abiertos). Se resaltan
          los dos umbrales fisicos reales: Xe en n_thr=4, Ar en n_thr=1.
Panel B : el integrando de kappa,  w(T) * dS/dEEE|_T  vs. energia de retroceso,
          en el umbral fisico de cada blanco -- muestra en que parte del
          espectro de retroceso vive la sensibilidad a EEE. Sombreado: T_thr
          (energia a la que <N_e> alcanza el borde del ROI).

Salida: datos/fig_kappa_eee.png
"""
import os
import re
import sys
import numpy as np
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kappa_eee import kappa_analytic, load_grid, EEE_NOM, EEE_SIGREL, NTHR_PHYS

DATOS = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
C = {"Xe": "#1f77b4", "Ar": "#E87722"}
NTHRS = (1, 2, 3, 4)

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
    "axes.linewidth": 0.9, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.fontsize": 8.6, "lines.linewidth": 1.7,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})


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


def t_thr(tag, n_thr):
    """energia de retroceso a la que <N_e creados> = n_thr (borde del ROI)."""
    T, w, nF, pF, lam = load_grid(tag)
    idx = np.where(lam >= n_thr)[0]
    return T[idx[0]] if len(idx) else np.nan


def main():
    fd = read_fd_kappa()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.3))

    # ---------- Panel A : kappa vs n_thr ----------
    a1.axhline(1.0, color="0.6", lw=1.0, ls="--")
    a1.text(1.02, 1.0, "normalización pura ($\\kappa=1$)", fontsize=7.6,
            color="0.45", va="bottom")
    for tag in ("Xe", "Ar"):
        ka = [kappa_analytic(tag, EEE_NOM[tag], n)[0] for n in NTHRS]
        a1.plot(NTHRS, ka, "-o", color=C[tag], ms=5, label=f"{tag} (analítico)")
        kf = [fd[tag].get(n, np.nan) for n in NTHRS]
        a1.plot(NTHRS, kf, "o", mfc="none", mec=C[tag], ms=10, mew=1.2,
                label=f"{tag} (dif. finita)")
        # umbral fisico
        nphys = NTHR_PHYS[tag]
        kphys = kappa_analytic(tag, EEE_NOM[tag], nphys)[0]
        a1.plot(nphys, kphys, marker="*", color=C[tag], ms=17, mec="black",
                mew=0.6, zorder=5)

    kx4 = kappa_analytic("Xe", EEE_NOM["Xe"], 4)[0]
    ka1 = kappa_analytic("Ar", EEE_NOM["Ar"], 1)[0]
    a1.annotate(rf"$n_{{\rm thr}}=4$ (Xe real): $\kappa={kx4:.2f}$" "\n"
                rf"$\kappa\cdot\sigma_{{\rm EEE}}/{{\rm EEE}}\approx{kx4*EEE_SIGREL['Xe']*100:.0f}\%$",
                xy=(4, kx4), xytext=(1.55, 3.05), fontsize=8.4, color=C["Xe"],
                arrowprops=dict(arrowstyle="->", color=C["Xe"], lw=0.8))
    a1.annotate(rf"$n_{{\rm thr}}=1$ (Ar real): $\kappa={ka1:.2f}$" "\n"
                r"(pequeño pero $\neq 0$)",
                xy=(1, ka1), xytext=(2.35, 0.18), fontsize=8.4, color=C["Ar"],
                ha="center",
                arrowprops=dict(arrowstyle="->", color=C["Ar"], lw=0.8))

    a1.set_xlabel(r"Umbral del ROI  $n_{\rm thr}$  ($N_e \geq n_{\rm thr}$)")
    a1.set_ylabel(r"$\kappa_{\rm EEE} = (\mathrm{EEE}/S)\;dS/d\mathrm{EEE}$")
    a1.set_title(r"Amplificación de la incertidumbre de EEE (comparación ideal)")
    a1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    a1.set_xlim(0.8, 4.4)
    a1.set_ylim(0, 3.8)
    a1.legend(loc="upper left", ncol=2, columnspacing=1.0, handletextpad=0.4)

    # ---------- Panel B : integrando w * dS/dEEE ----------
    for tag in ("Xe", "Ar"):
        nphys = NTHR_PHYS[tag]
        _, _, _, _, (T, w, nF, pF, p_eff, S_T, dS_T) = kappa_analytic(
            tag, EEE_NOM[tag], nphys)
        integ = w * dS_T
        # suaviza el escalonado de nint(<N_e>/p_F) para ver la envolvente
        integ = uniform_filter1d(np.maximum(integ, 0.0), size=9)
        integ = integ / integ.max()                 # normalizado al maximo
        a2.plot(T, integ, color=C[tag], lw=1.7,
                label=rf"{tag}, $n_{{\rm thr}}={nphys}$")
        tt = t_thr(tag, nphys)
        a2.axvline(tt, color=C[tag], lw=0.9, ls=":")
        a2.text(tt, 0.55, rf" $T_{{\rm thr}}^{{\rm {tag}}}$", color=C[tag],
                fontsize=8.0, rotation=90, va="bottom",
                transform=a2.get_xaxis_transform())

    a2.set_xlabel(r"Energía de retroceso nuclear  $T_{\rm nr}$  [keV]")
    a2.set_ylabel(r"$w(T)\;dS/d\mathrm{EEE}\,|_T$   (envolvente, normalizada al máx.)")
    a2.set_ylim(0, 1.15)
    a2.set_title(r"Dónde vive la sensibilidad a EEE en el espectro")
    a2.set_xlim(0.15, 1.6)
    a2.legend(loc="upper right")

    fig.tight_layout()
    out = f"{DATOS}/fig_kappa_eee.png"
    fig.savefig(out)
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
