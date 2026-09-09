#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Significancia de descubrimiento de Asimov, Z, complementaria a A_90.

A_90 responde "que amplitud excluiria esta exposicion". Z responde "cuan
lejos estamos de DETECTAR la senal SM":

    Z_A = sqrt( 2 [ (S+B) ln(1 + S/B) - S ] )        (B > 0)
    Z_A -> sqrt(2 S)  (limite Poisson B -> 0)

y, como S y B escalan linealmente con la exposicion, Z ~ sqrt(exposicion),
asi que la exposicion para NΣ es  mu_NΣ = mu_0 (NΣ / Z_0)^2.

Fuentes:
  Ar: datos/sensib_bkg_Ar.dat (S_tot, B_tot por escenario de fondo de 39Ar,
      exposicion fija 62 kg*dia, ref.[46]).
  Xe: sin fondo simulado -> B = 0, S = Ntot_ROI@192 de sensib_ideal_Xe.dat
      por umbral N_e >= k (exposicion base 192 kg*dia).

Salida: datos/fig_discovery_Z.png  +  datos/tabla_discovery_Z.tex
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
EXPO_XE = 192.0
EXPO_AR = 62.0
C_XE, C_AR = "#33546e", "#a86a43"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
    "axes.linewidth": 0.9, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "legend.fontsize": 8.4, "lines.linewidth": 1.7,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})


def Z_asimov(S, B):
    S = np.asarray(S, float); B = np.asarray(B, float)
    out = np.sqrt(2.0 * S)                      # limite B -> 0
    m = B > 1e-30
    out = np.where(m, np.sqrt(np.maximum(
        2.0 * ((S + B) * np.log1p(np.divide(S, B, where=m, out=np.zeros_like(S))) - S),
        0.0)), out)
    return out


def expo_for(mu0, Z0, Nsig):
    return mu0 * (Nsig / Z0) ** 2


# ---------------------------------------------------------------- Ar
ar = []
with open(f"{BASE}/sensib_bkg_Ar.dat", encoding="utf-8") as f:
    for ln in f:
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        ar.append((p[0], float(p[3]), float(p[4])))   # escenario, S_tot, B_tot

# ---------------------------------------------------------------- Xe ideal (B=0)
xe = {}
with open(f"{BASE}/sensib_ideal_Xe.dat", encoding="utf-8") as f:
    for ln in f:
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if p[1] == "Fnest":
            xe[int(p[0])] = float(p[4])                # Ntot_ROI@192

print("=" * 80)
print(" SIGNIFICANCIA DE DESCUBRIMIENTO (Z de Asimov)")
print("=" * 80)
rows_tex = []

print(f"\n Ar  (exposicion fija {EXPO_AR:.0f} kg*dia, ref.[46])")
print(f" {'escenario':22s} {'S_tot':>10s} {'B_tot':>12s} {'Z':>8s} "
      f"{'expo 3sig':>11s} {'expo 5sig':>11s}  [kg*dia]")
for name, S, B in ar:
    Z = float(Z_asimov(S, B))
    e3 = expo_for(EXPO_AR, Z, 3.0)
    e5 = expo_for(EXPO_AR, Z, 5.0)
    print(f" {name:22s} {S:10.1f} {B:12.3g} {Z:8.1f} {e3:11.2f} {e5:11.2f}")
    rows_tex.append(("Ar / " + name.replace("_", " "), S, B, Z, e3, e5))

print(f"\n Xe ideal  (B=0, exposicion base {EXPO_XE:.0f} kg*dia)")
print(f" {'umbral':22s} {'S=Ntot@192':>10s} {'B_tot':>12s} {'Z':>8s} "
      f"{'expo 3sig':>11s} {'expo 5sig':>11s}")
for k in (1, 2, 3, 4):
    S = xe[k]; Z = float(Z_asimov(S, 0.0))
    e3 = expo_for(EXPO_XE, Z, 3.0); e5 = expo_for(EXPO_XE, Z, 5.0)
    print(f" Xe N_e>={k:<15d} {S:10.1f} {'0':>12s} {Z:8.1f} {e3:11.1f} {e5:11.1f}")
    rows_tex.append((fr"Xe ideal $N_e\geq{k}$", S, 0.0, Z, e3, e5))

# ---------------------------------------------------------------- figura
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

# panel 1: Z por escenario (barras)
names = [r[0] for r in rows_tex]
Zs = [r[3] for r in rows_tex]
cols = [C_AR if n.startswith("Ar") else C_XE for n in names]
a1.barh(range(len(names)), Zs, color=cols, alpha=0.9, edgecolor=C_SM, lw=0.6)
a1.set_yticks(range(len(names)))
a1.set_yticklabels(names, fontsize=8)
a1.set_xscale("log")
for x, lab in ((3, r"$3\sigma$"), (5, r"$5\sigma$")):
    a1.axvline(x, color="0.4", lw=1.0, ls="--")
    a1.text(x, len(names) - 0.4, lab, fontsize=8, ha="center", va="bottom")
a1.set_xlabel(r"$Z$ de Asimov (a la exposición base)")
a1.set_title("Significancia de descubrimiento por escenario")
a1.invert_yaxis()

# panel 2: Z vs exposicion (Z = Z0 * sqrt(mu/mu0)), solo los casos con Z0 finito
# y no trivialmente enorme (los limitados por fondo o estadistica)
mu = np.logspace(0, 3.7, 200)
SHOW = {"Ar / UAr SEfloor ref46": (EXPO_AR, C_AR, "-"),
        r"Xe ideal $N_e\geq4$": (EXPO_XE, C_XE, "-"),
        r"Xe ideal $N_e\geq3$": (EXPO_XE, C_XE, ":")}
for name, S, B, Z0, _, _ in rows_tex:
    if name not in SHOW:
        continue
    mu0, col, ls = SHOW[name]
    a2.plot(mu, Z0 * np.sqrt(mu / mu0), color=col, ls=ls, label=name)
a2.axhline(3, color="0.4", lw=1.0, ls="--"); a2.axhline(5, color="0.4", lw=1.0, ls="--")
a2.text(1.3, 3, r"$3\sigma$", fontsize=8, va="bottom")
a2.text(1.3, 5, r"$5\sigma$", fontsize=8, va="bottom")
a2.set_xscale("log"); a2.set_yscale("log")
a2.set_xlabel("Exposición [kg·día]")
a2.set_ylabel(r"$Z$")
a2.set_title(r"$Z\propto\sqrt{\rm exposición}$: casos no triviales")
a2.legend(loc="lower right", fontsize=7.6)
a2.set_ylim(1, 60)

fig.suptitle("Descubrimiento: con Ar depletado es trivial ($Z\\!\\sim\\!160$); el "
             "único caso al borde es el piso de apilamiento SE ($Z\\!\\approx\\!4$) "
             "y Xe ideal con $N_e\\geq4$ ($Z\\!\\approx\\!5$)", y=1.02, fontsize=10.5)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_discovery_Z.png")

with open(f"{BASE}/tabla_discovery_Z.tex", "w", encoding="utf-8") as f:
    f.write("% python/discovery_Z.py\n\\begin{tabular}{@{}lccccc@{}}\n\\toprule\n")
    f.write("Caso & $S$ & $B$ & $Z$ & $\\mu_{3\\sigma}$ & $\\mu_{5\\sigma}$ "
            "[kg$\\cdot$día] \\\\\n\\midrule\n")
    for name, S, B, Z, e3, e5 in rows_tex:
        f.write(f"{name} & {S:.0f} & {B:.3g} & {Z:.1f} & {e3:.1f} & {e5:.1f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print(f"\n  {BASE}/fig_discovery_Z.png\n  {BASE}/tabla_discovery_Z.tex")
