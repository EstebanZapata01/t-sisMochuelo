#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lectura inversa de la curva de sensibilidad ideal: cuanta exposicion
necesita cada blanco (Xe, Ar) para alcanzar un A_90 de referencia, a
cada umbral N_e >= k. Es la pregunta complementaria a "que A_90 da esta
exposicion" (Tabla de sensib_ideal): aqui se fija el A_90 y se despeja
la exposicion.

Como la comparacion ideal es Asimov de conteo puro, A_90(mu) = 1 +
sqrt(2.706 / (mu * N_tot@192)); se invierte en forma cerrada:
    mu(A_ref) = 2.706 / (N_tot@192 * (A_ref - 1)^2)
    exposicion = 192 * mu   [kg*dia]
Se comprueba contra las 7 columnas A_90(x1..x335) de sensib_ideal_*.

Entradas: datos/sensib_ideal_{Xe,Ar}.dat  (fila F_mode = Fnest)
Salidas : stdout + datos/tabla_expo_A90.tex
"""
import numpy as np

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
EXPO_BASE = 192.0
DCHI2 = 2.706
A_REF = [1.10, 1.05, 1.02]
MULT = np.array([1, 2, 5, 10, 50, 100, 335], dtype=float)


def load_fnest(path):
    """-> dict  NE_LO -> (Ntot@192, A_90[x1..x335])."""
    out = {}
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            if p[1] != "Fnest":
                continue
            ne = int(p[0])
            ntot = float(p[4])
            a90 = np.array([float(x) for x in p[5:12]])
            out[ne] = (ntot, a90)
    return out


xe = load_fnest(f"{BASE}/sensib_ideal_Xe.dat")
ar = load_fnest(f"{BASE}/sensib_ideal_Ar.dat")


def expo_for(ntot, a_ref):
    mu = DCHI2 / (ntot * (a_ref - 1.0) ** 2)
    return EXPO_BASE * mu


rows = []
print("=" * 88)
print(" Exposición [kg·día] para alcanzar A_90 de referencia (comparación ideal)")
print("=" * 88)
hdr = f"{'N_e>=':>6s} " + "".join(f"{'Xe @'+f'{a:.2f}':>14s}{'Ar @'+f'{a:.2f}':>14s}" for a in A_REF)
print(hdr)
print("-" * 88)
for ne in (1, 2, 3, 4):
    ntx, _ = xe[ne]
    nta, _ = ar[ne]
    cells = []
    for a in A_REF:
        ex, ea = expo_for(ntx, a), expo_for(nta, a)
        cells.append((ex, ea))
    print(f"{ne:>6d} " + "".join(f"{ex:14.0f}{ea:14.0f}" for ex, ea in cells))
    rows.append((ne, cells))
print("-" * 88)

# chequeo: la formula cerrada debe reproducir las columnas del .dat
ne = 4
ntx, a90x = xe[ne]
pred = 1.0 + np.sqrt(DCHI2 / (MULT * ntx))
print(f" chequeo N_e>=4 Xe:  A_90(dat)  = {np.array2string(a90x, precision=3)}")
print(f"                     A_90(form) = {np.array2string(pred, precision=3)}")

with open(f"{BASE}/tabla_expo_A90.tex", "w", encoding="utf-8") as f:
    f.write("% generado por python/expo_to_A90.py  (exposicion en kg*dia)\n")
    f.write("\\begin{tabular}{@{}c" + "cc" * len(A_REF) + "@{}}\n\\toprule\n")
    f.write("$N_e\\ge$ & " + " & ".join(
        f"\\multicolumn{{2}}{{c}}{{$A_{{90}}={a:.2f}$}}" for a in A_REF) + " \\\\\n")
    f.write(" & " + " & ".join("Xe & Ar" for _ in A_REF) + " \\\\\n\\midrule\n")
    for ne, cells in rows:
        f.write(f"{ne} & " + " & ".join(
            f"{ex:.0f} & {ea:.0f}" for ex, ea in cells) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print(f"\n  {BASE}/tabla_expo_A90.tex")
