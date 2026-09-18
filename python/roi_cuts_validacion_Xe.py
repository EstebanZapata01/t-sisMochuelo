#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe, ROI (N_e=4-7): compara la simulacion contra las dos figuras del paper
que muestran el espectro en N_e -- Fig. 3 ("N_e extraidos", sin cortes de
seleccion) y Fig. 6 ("CEvNS signal before/after cuts") -- para responder si
son la misma curva y si la simulacion reproduce alguna de las dos.

Panel A (antes de cortes): sim = tasa_ion_extraidos (misma cantidad que
  entra en chi2.f90 antes de eff_ROI). Sigue la FORMA de la Fig. 3 (razones
  entre bins consecutivos ~6-7 en ambas) con una normalizacion ~1.8-2.1x
  mayor (espectro hibrido mas duro que SM2018, ya documentado). NO sigue a
  la Fig. 6 "before cuts": esa curva es mucho mas plana (razon 4->5 = 1.2,
  no ~7) y ~2-4x mas baja que la simulacion en cada bin. Fig. 3 y Fig. 6
  "before cuts" NO son la misma curva dentro del propio paper.
Panel B (despues de cortes): sim*eff_ROI vs Fig. 6 "after cuts". eff_ROI se
  digitalizo como razon Fig6_despues/Fig6_antes, asi que aplicarla sobre una
  forma tipo Fig. 3 (no tipo Fig. 6) no la reconcilia bin a bin: N_e=4 sale
  sobre-predicho ~3.8x, N_e=5-7 sub-predichos ~0.6-0.7x -- ya documentado en
  metodologia.tex. Lo que falta es la aceptancia del detector (eficiencia de
  trigger/borde cerca de N_e=4, duracion de cluster, radio^2 reconstruido),
  que ni la Fig. 3 ni la simulacion tienen y que solo el Monte Carlo de
  RED-100 podria proveer.

Entrada: datos/validacion_fig3_fig6_Xe.dat (volcado de red100PE.f90; valores
         Fig. 3 y Fig. 6 son datos digitalizados por el usuario del paper,
         arXiv:2411.18641).
Salida : datos/fig_roi_cuts_Xe.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, GRIS, NEGRO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

d = np.loadtxt(f"{BASE}/validacion_fig3_fig6_Xe.dat", comments="#")
Ne, sim_a, sim_d, fig3, f6_a, f6_d = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
x = np.arange(len(Ne))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.6, 4.6))

w = 0.26
axA.bar(x - w, sim_a, width=w, color=C_XE, label="Simulación (NEST)")
axA.bar(x, fig3, width=w, color=GRIS, label="RED-100, antes de reconstrucción")
axA.bar(x + w, f6_a, width=w, color=NEGRO, label="RED-100, tras reconstrucción")
axA.set_yscale("log")
axA.set_xticks(x); axA.set_xticklabels([f"{int(n)}" for n in Ne])
axA.set_xlabel(r"$N_e$")
axA.set_ylabel(r"eventos / (kg$\cdot$día)")
axA.set_title("Antes de cortes")
axA.legend(loc="upper right", fontsize=7.6)

w2 = 0.32
axB.bar(x - w2/2, sim_d, width=w2, color=C_XE, label="Simulación $\\times\\,$eff$_{\\rm ROI}$")
axB.bar(x + w2/2, f6_d, width=w2, color=NEGRO, label="RED-100, tras cortes")
axB.set_yscale("log")
axB.set_xticks(x); axB.set_xticklabels([f"{int(n)}" for n in Ne])
axB.set_xlabel(r"$N_e$")
axB.set_title("Después de cortes")
axB.legend(loc="upper right", fontsize=7.6)

fig.tight_layout()
out = f"{BASE}/fig_roi_cuts_Xe.png"
fig.savefig(out)

print("=" * 70)
print(" Xe, ROI: simulacion vs Fig. 3 vs Fig. 6 (antes/despues de cortes)")
print("=" * 70)
for i in range(len(Ne)):
    print(f"  Ne={int(Ne[i])}: sim/Fig3={sim_a[i]/fig3[i]:.2f}  "
          f"sim/Fig6antes={sim_a[i]/f6_a[i]:.2f}  "
          f"(sim*eff)/Fig6despues={sim_d[i]/f6_d[i]:.2f}")
print(f"\n  {out}")
