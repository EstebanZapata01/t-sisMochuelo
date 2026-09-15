#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tabla-resumen de los tres blancos CEvNS considerados (Xe, Ge, Ar):
masa, (Z, N), N/Z, umbral de retroceso nuclear T_thr [keV_nr], Q_W, Q_W^2,
radio del anillo ciego rho_eps y direccion ciega theta.

T_thr se define como la energia de retroceso a la que el numero medio de
electrones de ionizacion CREADOS, <N_e> = T * Qy(T), alcanza el borde
inferior del ROI de ese analisis:
    Xe : N_e >= 4   (RED-100, corte por fondo de electron unico)
    Ar : N_e >= 1   (RED-100 §VII, "below five ionization electrons")
    Ge : E_er >= 160 eV_ee  (CONUS+); se convierte a T_nr via quenching
         de Lindhard Q(T), T_thr donde Q(T)*T = 0.160 keV.

Entradas: datos/nest_218V_dense.txt, datos/nest_Ar_218V_dense.txt
Salidas : stdout + datos/tabla_umbral.tex
"""
import numpy as np

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

# ------------------------------------------------------------------ blancos
# Q_W del SM con sin^2(theta_W) = 0.23857  (Q_W = -N/2 + 0.02286 Z)
S2W = 0.23857
def QW(Z, N):
    return -N / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * Z

TARGETS = {
    "Xe-131": dict(Z=54, N=77, A=131, roi_ne=4),
    "Ge-73":  dict(Z=32, N=41, A=73,  roi_ne=None),   # umbral por energia
    "Ar-40":  dict(Z=18, N=22, A=40,  roi_ne=1),
}

# ------------------------------------------------------------------ T_thr NEST
def tthr_from_ne(table_path, ne_target):
    T, Qy = np.loadtxt(table_path, comments="#", usecols=(0, 1), unpack=True)
    ne = T * Qy                       # <N_e> creados
    # primer T donde <N_e> cruza ne_target
    idx = np.argmax(ne >= ne_target)
    if idx == 0:
        return np.nan
    # interpolacion lineal en el cruce
    t0, t1 = T[idx - 1], T[idx]
    n0, n1 = ne[idx - 1], ne[idx]
    return t0 + (ne_target - n0) * (t1 - t0) / (n1 - n0)

tthr_xe = tthr_from_ne(f"{BASE}/nest_218V_dense.txt", 4)
tthr_ar = tthr_from_ne(f"{BASE}/nest_Ar_218V_dense.txt", 1)

# ------------------------------------------------------------------ T_thr Ge (Lindhard)
def Q_lindhard(T_keV, Z=32):
    eps = 11.5 * Z ** (-7.0 / 3.0) * T_keV
    g = 3.0 * eps ** 0.15 + 0.7 * eps ** 0.6 + eps
    k = 0.162
    return k * g / (1.0 + k * g)

def tthr_ge(Eer_keV=0.160):
    T = np.linspace(1e-3, 5.0, 200000)
    Eer = Q_lindhard(T) * T
    idx = np.argmax(Eer >= Eer_keV)
    return T[idx]

tthr_ge = tthr_ge()

TTHR = {"Xe-131": tthr_xe, "Ge-73": tthr_ge, "Ar-40": tthr_ar}

# ------------------------------------------------------------------ tabla
rows = []
print("=" * 92)
print(f"{'blanco':8s} {'Z':>3s} {'N':>3s} {'N/Z':>6s} {'T_thr[keV_nr]':>13s} "
      f"{'Q_W':>9s} {'Q_W^2':>9s} {'rho_eps':>8s} {'theta[deg]':>11s}")
print("-" * 92)
for name, d in TARGETS.items():
    Z, N = d["Z"], d["N"]
    qw = QW(Z, N)
    rho = abs(qw) / (Z + 2 * N)
    theta = np.degrees(np.arctan2(Z + 2 * N, 2 * Z + N))
    tt = TTHR[name]
    print(f"{name:8s} {Z:3d} {N:3d} {N/Z:6.3f} {tt:13.3f} "
          f"{qw:9.3f} {qw**2:9.1f} {rho:8.4f} {theta:11.3f}")
    rows.append((name, Z, N, N / Z, tt, qw, qw**2, rho, theta))
print("-" * 92)
th = {r[0]: r[8] for r in rows}
print(f"  angulo ciego Xe-Ar = {abs(th['Xe-131'] - th['Ar-40']):.2f} deg")
print(f"  angulo ciego O-Xe (maximo posible, O-16 Z=N=8) = "
      f"{abs(np.degrees(np.arctan2(24,24)) - th['Xe-131']):.2f} deg")

# ------------------------------------------------------------------ .tex
with open(f"{BASE}/tabla_umbral.tex", "w", encoding="utf-8") as f:
    f.write("% generado por python/threshold_table.py\n")
    f.write("\\begin{tabular}{@{}lccccccc@{}}\n\\toprule\n")
    f.write("Blanco & $Z$ & $N$ & $N/Z$ & $T_{\\rm thr}$ [keV$_{nr}$] & "
            "$\\QW$ & $\\QW^{2}$ & $\\theta$ [$^\\circ$] \\\\\n\\midrule\n")
    labels = {"Xe-131": "$^{131}$Xe (RED-100)",
              "Ge-73": "$^{73}$Ge (CONUS+)",
              "Ar-40": "$^{40}$Ar (RED-100)"}
    for name, Z, N, nz, tt, qw, qw2, rho, theta in rows:
        f.write(f"{labels[name]} & {Z} & {N} & {nz:.3f} & {tt:.2f} & "
                f"{qw:.2f} & {qw2:.1f} & {theta:.2f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print(f"\n  {BASE}/tabla_umbral.tex")
