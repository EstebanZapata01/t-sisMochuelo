#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validacion de Ar contra la ref.[46] (Physics 5, 492 (2023)), SIN comparar
contra datos reales de Xe y SIN barrer la exposicion.

El observable fisico es el mismo que usa la propia ref.[46]: S/sqrt(B)
(relacion senal-fondo). Se calcula solo a partir de la prediccion SM (S_k,
CEvNS) y del fondo simulado (B_k, de actividades de 39Ar PUBLICADAS) en la
UNICA exposicion real y justificada para Ar: 62 kg (masa del "proximo
montaje" de ref.[46]) x 1 dia = 62 kg*dia -- el mismo numero que ref.[46]
declara para su S/sqrt(B)~4. No hay ningun plan publicado de correr Ar a
otra exposicion, asi que no se barre nada.

El "parametro que se mueve" es el escenario de FONDO (un error/incertidumbre
discreta y justificada), no la exposicion -- cuatro corridas de reemplazo
("Sec. brecha ideal->real" de metodologia.tex), calculadas por
chi2_bkg_nest.f90 y escritas en datos/sensib_bkg_Ar.dat:
    sin_fondo          -- techo intrinseco (sin ningun fondo)
    Ar39_UAr           -- 39Ar depletado (DarkSide-50, 7.3e-4 Bq/kg)
    Ar39_atmosferico   -- 39Ar sin depletar (ref.[46], 1 Bq/kg)
    UAr_SEfloor_ref46  -- UAr + el piso de apilamiento de electron unico
                          CALIBRADO para reproducir el S/sqrt(B)~4 que
                          ref.[46] declara (ese fondo lo deja sin resolver;
                          se infiere su tamano implicito del propio numero
                          publicado, no se supone nada mas)

Salida: datos/fig_ar_fondo_validacion.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
SB_REF46 = 4.0

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9, "axes.grid": True,
    "grid.color": "0.82", "grid.linewidth": 0.5,
    "xtick.direction": "in", "ytick.direction": "in",
    "ytick.right": True,
    "legend.frameon": False, "legend.fontsize": 8.5,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})

LABELS = {
    "sin_fondo": "Sin fondo\n(techo intrínseco)",
    "Ar39_UAr": "$^{39}$Ar UAr\n(depletado)",
    "Ar39_atmosferico": "$^{39}$Ar\natmosférico",
    "UAr_SEfloor_ref46": "UAr + piso SE\n(calibrado a ref.[46])",
}
ORDER = ["sin_fondo", "Ar39_UAr", "Ar39_atmosferico", "UAr_SEfloor_ref46"]

rows = {}
expo = None
with open(f"{BASE}/sensib_bkg_Ar.dat", encoding="utf-8") as f:
    for ln in f:
        if ln.lstrip().startswith("#") or not ln.strip():
            continue
        p = ln.split()
        esc, e, a90, s_tot, b_tot, sb = p[0], float(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[5])
        rows[esc] = dict(expo=e, a90=a90, S=s_tot, B=b_tot, sb=sb)
        expo = e

print("=" * 78)
print(f" VALIDACION DE Ar CONTRA ref.[46]  (exposicion UNICA y real: {expo:.0f} kg*dia)")
print("=" * 78)
for esc in ORDER:
    r = rows[esc]
    print(f"   {esc:20s}  A_90={r['a90']:7.4f}   S/sqrt(B)={r['sb']:9.3f}")
print(f" ref.[46] declara S/sqrt(B) ~ {SB_REF46:.1f} a la misma exposicion -- el escenario")
print(" 'UAr_SEfloor_ref46' se CALIBRA para reproducirlo exactamente (no es una")
print(" suposicion nueva: es el tamano implicito del fondo que ref.[46] deja sin")
print(" resolver, inferido de su propio numero publicado).")

# "sin_fondo" no tiene S/sqrt(B) definido (B=0, techo intrinseco): se anota
# aparte, no como barra en un eje log donde una barra de altura 0 es invisible.
ORDER_SB = ["Ar39_UAr", "Ar39_atmosferico", "UAr_SEfloor_ref46"]
sb_vals = [rows[e]["sb"] for e in ORDER_SB]
a90_vals = [rows[e]["a90"] for e in ORDER_SB]
labels = [LABELS[e] for e in ORDER_SB]
x = np.arange(len(ORDER_SB))

fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.bar(x, sb_vals, color="0.75", edgecolor="black", width=0.55, zorder=2)
ax.axhline(SB_REF46, color="black", lw=1.3, ls="--", zorder=3)
ax.text(0.05, SB_REF46 * 1.4, r"ref.[46] declara $S/\sqrt{B}\approx 4$",
        transform=ax.get_yaxis_transform(), fontsize=8.5, ha="left", va="bottom")
for xi, sb, a90 in zip(x, sb_vals, a90_vals):
    ax.annotate(f"$A_{{90}}={a90:.2f}$", (xi, sb), textcoords="offset points",
                xytext=(0, 6), ha="center", fontsize=8.3)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.7)
ax.set_yscale("log")
ax.set_ylabel(r"$S/\sqrt{B}$  (observable de ref.[46])")
ax.set_title(f"Ar: validación contra ref.[46] a su propia exposición real "
             f"({expo:.0f} kg$\\cdot$día)")
ax.set_ylim(1, 1e4)
ax.text(0.985, 0.60, f"(sin fondo: techo intrínseco, $A_{{90}}={rows['sin_fondo']['a90']:.2f}$,\n"
        r"$S/\sqrt{B}$ no definido — no graficado)",
        transform=ax.transAxes, fontsize=7.6, ha="right", va="bottom", color="0.35")
fig.tight_layout()
fig.savefig(f"{BASE}/fig_ar_fondo_validacion.png")
print(f"\n  {BASE}/fig_ar_fondo_validacion.png")
