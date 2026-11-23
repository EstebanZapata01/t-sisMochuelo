#!/usr/bin/env python3
"""
Tabla de NEST para el xenon (retroceso nuclear, 218 V/cm) con nestpy.
Salida: datos/nest_218V_dense.txt, que lee Tnr_to_e.f90 (carpeta Xe):
    T_nr[keV]   Qy[e-/keV]   F = Var(N_e)/<N_e>

Qy es el rendimiento medio de carga y F se estima muestreando GetQuanta
N_SAMP veces por energia. Condiciones de RED-100 (arXiv:2411.18641):
169 K, 218 V/cm, densidad 2.96 g/cm3.
"""
from pathlib import Path
import nestpy
import numpy as np

# Semilla fija para que la tabla sea reproducible
nestpy.RandomGen.rndm().set_seed(20260910)
nestpy.RandomGen.rndm().lock_seed()

NR = nestpy.interactions.NR
DENSITY   = 2.96      # g/cm3, LXe a ~169 K
DRIFT_V   = 218.0     # V/cm
N_SAMP    = 20000     # muestras de GetQuanta por energia para estimar F
T_MIN, T_MAX, N_T = 0.2, 12.0, 4000   # el retroceso maximo de CEvNS de reactor es ~1 keV

print(f"nestpy / NEST: {nestpy.__nest_version__}")

detector = nestpy.detectors.DetectorExample_XENON10()
nc = nestpy.NESTcalc(detector)

# Parametros por defecto de NEST v2.4.0 para NR en Xe
nr_params = [11.0, 1.1, 0.048, -0.0533, 12.6, 0.3, 2.0, 0.3, 2.0, 0.5, 1.0, 1.0]
er_params = [-1.0] * 10

T_vals  = np.linspace(T_MIN, T_MAX, N_T)
Qy_vals = np.zeros(N_T)
F_vals  = np.ones(N_T)          # F = 1 (Poisson) como valor de respaldo

print(f"Calculando Qy y F ({N_T} energias x {N_SAMP} muestras)...")
for i, T in enumerate(T_vals):
    try:
        y = nc.GetYields(interaction=NR, energy=float(T), density=DENSITY,
                         drift_field=DRIFT_V, A=131, Z=54,
                         nr_parameters=nr_params, er_parameters=er_params)
        Qy_vals[i] = y.ElectronYield / T

        ne = np.fromiter((nc.GetQuanta(y, DENSITY).electrons for _ in range(N_SAMP)),
                         dtype=float, count=N_SAMP)
        m = ne.mean()
        if m > 1e-6:
            F_vals[i] = ne.var() / m
    except Exception as e:
        print(f"  Error en T={T:.3f}: {e}")
        Qy_vals[i] = 0.0
        F_vals[i]  = 1.0

output = str(Path(__file__).resolve().parents[1] / "datos" / "nest_218V_dense.txt")
with open(output, "w") as f:
    f.write("# Tabla NEST para RED-100\n")
    f.write(f"# Xe NR, rho={DENSITY} g/cm3, E_drift={DRIFT_V} V/cm, NEST {nestpy.__nest_version__}\n")
    f.write(f"# F estimado con {N_SAMP} muestras de GetQuanta por energia\n")
    f.write("# T_nr[keV]   Qy[e-/keV]   F=Var(Ne)/<Ne>\n")
    f.write(f"# N_puntos = {N_T}\n")
    for T, Qy, F in zip(T_vals, Qy_vals, F_vals):
        f.write(f"{T:.6f}   {Qy:.6f}   {F:.6f}\n")

print(f"\nArchivo '{output}' creado. Verificacion:")
for T_test in [0.21, 0.5, 1.0, 3.0, 5.0]:
    idx = int(np.argmin(np.abs(T_vals - T_test)))
    print(f"  T={T_vals[idx]:6.3f} keV -> Qy={Qy_vals[idx]:.4f} e-/keV   F={F_vals[idx]:.4f}")
