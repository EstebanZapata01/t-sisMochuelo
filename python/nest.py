#!/usr/bin/env python3
"""
Genera la tabla NEST para RED-100 (Xe NR, 218 V/cm) con nestpy 2.4.5.

Tesis: metodologia.tex Sec. 3.2 (modelo binomial de Fano). Salida ->
       datos/nest_218V_dense.txt, leida por Tnr_to_e.f90 (carpeta Xe).


Salida: 3 columnas
    T_nr[keV]   Qy[e-/keV]   F = Var(N_e)/<N_e>   (factor tipo Fano)

- Qy  : rendimiento de carga medio  (y.ElectronYield / T)
- F   : ancho relativo de la distribucion REAL de N_e de NEST, obtenido
        muestreando GetQuanta N_SAMP veces por energia. Para NR en LXe
        F ~ 0.43-0.47 (sub-Poissoniano; Poisson seria F = 1).
        El Fortran usa (Qy, F) para modelar N_e creados como
        Binomial(n_F, p_F) con p_F = 1 - F, n_F = nint(<N_e>/p_F),
        y los extraidos como Binomial(n_F, p_F * EEE).

Parametros del paper arXiv:2411.18641 (Sec. II.B): NEST v2.4.0, 169 K,
1.29 bar, 218 V/cm.  rho(LXe, 169 K) ~ 2.96 g/cm3.
Los nr_parameters son los defaults de NEST v2.4.0 para NR en Xe.
"""
import nestpy
import numpy as np

NR = nestpy.interactions.NR
DENSITY   = 2.96      # g/cm3, LXe a ~169 K
DRIFT_V   = 218.0     # V/cm
N_SAMP    = 20000     # muestras de GetQuanta por energia para estimar F
T_MIN, T_MAX, N_T = 0.2, 12.0, 4000   # el retroceso maximo de CEvNS de reactor es ~1 keV

print(f"nestpy / NEST: {nestpy.__nest_version__}")

detector = nestpy.detectors.DetectorExample_XENON10()
nc = nestpy.NESTcalc(detector)

# Defaults de NEST v2.4.0 para NR en Xe (sacados del mensaje de C++)
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

output = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/nest_218V_dense.txt"
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
