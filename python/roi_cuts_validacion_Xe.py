#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Xe, ROI (N_e=4-7): compara la simulacion contra las dos figuras del paper
que muestran el espectro en N_e -- Fig. 3 ("N_e extraidos", sin cortes de
seleccion) y Fig. 6 ("CEvNS signal before/after cuts") -- para responder si
son la misma curva y si la simulacion reproduce alguna de las dos.

Etapas mostradas (2 series de "Simulacion", 3 de "RED-100"):
  teorico puro (Fig. 3 / sim = tasa_ion_extraidos): la simulacion SI sigue la
    forma de esta etapa (~1.8-2.1x mas alta, espectro hibrido mas duro que
    SM2018, ya documentado).
  tras reconstruccion, antes de cortes finales (Fig. 6 "before cuts"): solo
    RED-100 -- la simulacion no tiene esta etapa (necesita reconstruccion de
    posicion/duracion/energia, no reproducible solo con el paper).
  final, tras cortes (Fig. 6 "after cuts" vs sim*eff_ROI): N_e=4 sale
    sobre-predicho ~3.8x, N_e=5-7 sub-predichos ~0.6-0.7x -- eff_ROI se
    digitalizo de la forma de Fig. 6, no de Fig. 3, asi que aplicarla sobre
    una forma tipo Fig. 3 no reconcilia bin a bin. Falta la aceptancia del
    detector (eficiencia de trigger/borde cerca de N_e=4), que ni la Fig. 3
    ni la simulacion tienen.

Entrada: datos/validacion_fig3_fig6_Xe.dat (volcado de red100PE.f90; valores
         Fig. 3 y Fig. 6 son datos digitalizados por el usuario del paper,
         arXiv:2411.18641).
Salida : datos/fig_roi_cuts_Xe.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, NEGRO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

d = np.loadtxt(f"{BASE}/validacion_fig3_fig6_Xe.dat", comments="#")
Ne, sim_a, sim_d, fig3, f6_a, f6_d = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]

fig, ax = plt.subplots(figsize=(7.0, 5.2))

series = [
    (fig3,  NEGRO, "o", "RED-100, teórico"),
    (f6_a,  NEGRO, "^", "RED-100, tras reconstrucción"),
    (f6_d,  NEGRO, "s", "RED-100, tras cortes"),
    (sim_a, C_XE,  "o", "Simulación, teórico"),
    (sim_d, C_XE,  "s", "Simulación $\\times\\,$eff$_{\\rm ROI}$"),
]
for y, col, mk, lab in series:
    ax.plot(Ne, y, color=col, alpha=0.3, lw=1.2, zorder=2)
    ax.scatter(Ne, y, color=col, marker=mk, s=42, label=lab, zorder=3)

ax.set_yscale("log")
ax.set_xticks(Ne); ax.set_xticklabels([f"{int(n)}" for n in Ne])
ax.set_xlim(3.5, 7.5)
ax.set_xlabel(r"$N_e$")
ax.set_ylabel(r"eventos / (kg$\cdot$día)")
ax.legend(loc="upper right", fontsize=7.6)

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
