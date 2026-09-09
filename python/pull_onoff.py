#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Histograma de "pulls" (residuos normalizados) del ajuste ON-OFF real de Xe.

Para cada bin i:  pull_i = (dN_i - A_best * R_SM_i) / sigma_i.
Con A_best = 0 (senal CEvNS no requerida por los datos) esto es simplemente
dN_i / sigma_i. Si el residuo ON-OFF es ruido puro compatible con cero, los
15 pulls deben distribuirse como una N(0,1): media ~ 0, desviacion ~ 1.
Es la figura estandar en fisica de particulas para argumentar VISUALMENTE
"esto es ruido, no senal", mas contundente que reportar solo chi2/ndof.

Entrada: datos/chi2_ON_OFF_banda.dat
    columnas: PE_center  R_SM  A90*R_SM  -A90*R_SM  delta_ON_OFF  sigma_stat
    comentarios: "# A_best = ...", "# A_90 = ..."

Salida: datos/fig_pull_onoff.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()
from scipy import stats

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.fontsize": 8.5,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})

A_best = 0.0
rows = []
with open(f"{BASE}/chi2_ON_OFF_banda.dat", encoding="utf-8") as f:
    for ln in f:
        s = ln.strip()
        if s.startswith("#"):
            if "A_best" in s:
                A_best = float(s.split("=")[1])
            continue
        if not s:
            continue
        rows.append([float(x) for x in s.split()])
a = np.array(rows)
pe, R_SM, dN, sigma = a[:, 0], a[:, 1], a[:, 4], a[:, 5]

pull = (dN - A_best * R_SM) / sigma
mu, sd = pull.mean(), pull.std(ddof=1)
ks_stat, ks_p = stats.kstest(pull, "norm")
chi2_ndof = np.sum(pull ** 2) / len(pull)

print("=" * 70)
print(f" PULLS del ajuste ON-OFF (Xe real, A_best = {A_best:.3f})")
print("=" * 70)
print(f"  n bins        = {len(pull)}")
print(f"  media         = {mu:+.3f}   (esperado 0)")
print(f"  desv. tipica  = {sd:.3f}   (esperado 1)")
print(f"  sum pull^2/n  = {chi2_ndof:.3f}   (= chi2/ndof)")
print(f"  K-S vs N(0,1) = D {ks_stat:.3f},  p {ks_p:.2f}")

fig, (axh, axs) = plt.subplots(
    1, 2, figsize=(9.4, 4.0), gridspec_kw={"width_ratios": [1.0, 1.25]})

# --- panel izq: histograma vs N(0,1)
bins = np.linspace(-3, 3, 13)
axh.hist(pull, bins=bins, density=True, color="0.78", edgecolor="black",
         lw=0.8, label=f"pulls ({len(pull)} bins)")
xx = np.linspace(-3.4, 3.4, 300)
axh.plot(xx, stats.norm.pdf(xx), color="#a86a43", lw=1.8, label=r"$\mathcal{N}(0,1)$")
axh.axvline(0.0, color="0.5", lw=0.7, ls=":")
axh.set_xlabel(r"pull $=(\Delta N_i - A_{\rm best}R_i)/\sigma_i$")
axh.set_ylabel("densidad")
axh.set_title("Distribución de residuos normalizados")
axh.legend(loc="upper left")
axh.text(0.97, 0.95,
         f"media $={mu:+.2f}$\n"
         f"$\\sigma\\;\\,={sd:.2f}$\n"
         f"K-S $p={ks_p:.2f}$",
         transform=axh.transAxes, ha="right", va="top", fontsize=8.3,
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.7))

# --- panel der: pull por bin
axs.axhline(0.0, color="0.5", lw=0.8, ls="-")
for lv in (1, 2):
    axs.axhspan(-lv, lv, color="0.90" if lv == 2 else "0.80", zorder=0)
axs.plot(pe, pull, "o", color="black", ms=5, zorder=3)
axs.set_xlabel("Energía corregida [PE]")
axs.set_ylabel("pull por bin")
axs.set_title("Residuo normalizado vs. energía")
axs.set_ylim(-3.2, 3.2)
axs.text(0.02, 0.96, r"bandas: $\pm1\sigma$, $\pm2\sigma$", transform=axs.transAxes,
         ha="left", va="top", fontsize=8.0, color="0.35")

fig.suptitle("Xe (RED-100): el residuo ON$-$OFF real es compatible con ruido",
             y=1.02, fontsize=12)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_pull_onoff.png")
print(f"\n  {BASE}/fig_pull_onoff.png")
