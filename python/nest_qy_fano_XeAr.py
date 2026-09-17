#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rendimiento de carga Q_y(T) y factor de Fano F(T)=Var(N_e)/<N_e> de NEST (Xe)
y LArNEST (Ar), superpuestos, hasta T~2 keV (la ROI CEvNS de ambos blancos).
Muestra por que la fluctuacion de ionizacion del argon es mas benigna
(F_Ar ~ 0.11-0.15 vs F_Xe ~ 0.43-0.56).

Entradas: datos/nest_218V_dense.txt, datos/nest_Ar_218V_dense.txt
          (python/nest.py, python/nest_Ar.py; 3 columnas T_nr Qy F)
Salida  : datos/fig_nest_qy_fano.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
T_MAX = 2.0   # keV

dXe = np.loadtxt(f"{BASE}/nest_218V_dense.txt", comments="#")
dAr = np.loadtxt(f"{BASE}/nest_Ar_218V_dense.txt", comments="#")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 4.4))

for d, col, tag in ((dXe, C_XE, "Xe"), (dAr, C_AR, "Ar")):
    sel = d[:, 0] <= T_MAX
    a1.plot(d[sel, 0], d[sel, 1], color=col, label=tag)
    a2.plot(d[sel, 0], d[sel, 2], color=col, label=tag)

a1.set_xlim(0, T_MAX)
a1.set_ylim(0, None)
a1.set_xlabel(r"$T_{\rm nr}$  [keV]")
a1.set_ylabel(r"$Q_y(T_{\rm nr})$  [e$^-$/keV]")
a1.legend(loc="lower right")

a2.set_xlim(0, T_MAX)
a2.set_ylim(0, 0.7)
a2.set_xlabel(r"$T_{\rm nr}$  [keV]")
a2.set_ylabel(r"$F(T_{\rm nr}) = \mathrm{Var}(N_e)/\langle N_e\rangle$")
a2.legend(loc="center right")

fig.tight_layout()
out = f"{BASE}/fig_nest_qy_fano.png"
fig.savefig(out)

print("=" * 60)
for d, tag in ((dXe, "Xe"), (dAr, "Ar")):
    sel = d[:, 0] <= T_MAX
    print(f"  {tag}: F en [0,{T_MAX}] keV -> min={d[sel,2].min():.3f}  "
          f"max={d[sel,2].max():.3f}")
print(f"\n  {out}")
