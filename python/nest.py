#!/usr/bin/env python3
"""
Script para generar Qy con nestpy 2.4.5 usando DetectorExample_XENON10
y pasando todos los argumentos que pide GetYields.
"""
import nestpy
import numpy as np

print(f"nestpy versión: {nestpy.__nest_version__}")
print("Usando detector de ejemplo y GetYields con todos los parámetros.\n")

# Crear detector (el primer disponible)
detector = nestpy.detectors.DetectorExample_XENON10()
nc = nestpy.NESTcalc(detector)

# Valores por defecto de nr_parameters y er_parameters (sacados del error de C++)
nr_params = [11.0, 1.1, 0.048, -0.0533, 12.6, 0.3, 2.0, 0.3, 2.0, 0.5, 1.0, 1.0]
er_params = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

T_vals = np.linspace(0.2, 15, 10000)
Qy_vals = np.zeros(len(T_vals))

print("Calculando (puede tardar unos segundos)...")
for i, T in enumerate(T_vals):
    try:
        # CORRECCIÓN: Se usa nestpy.interactions.NR 
        y = nc.GetYields(interaction = nestpy.interactions.NR,
                         energy = T,
                         density = 2.9,
                         drift_field = 218.0,
                         A = 131,
                         Z = 54,
                         nr_parameters = nr_params,
                         er_parameters = er_params)
        
        # y.ElectronYield da el número total de electrones.
        # Para obtener Q_y en (e-/keV), lo dividimos por la energía T.
        Qy_vals[i] = y.ElectronYield / T
        
    except Exception as e:
        print(f"Error en T={T:.3f}: {e}")
        Qy_vals[i] = 0.0

# Guardar archivo en la ruta absoluta especificada
output = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/nest_218V_dense.txt"

with open(output, "w") as f:
    f.write("# Q_y tabla NEST para RED-100\n")
    f.write("# Xe NR, rho=2.9 g/cm3, E=218 V/cm\n")
    f.write("# T_nr[keV]   Qy[e-/keV]\n")
    f.write(f"# N_puntos = {len(T_vals)}\n")
    for T, Qy in zip(T_vals, Qy_vals):
        f.write(f"{T:.6f}   {Qy:.6f}\n")

print(f"\nArchivo '{output}' creado. Verificación:")
for T_test in [0.1, 0.2, 0.21, 3.0, 5]:
    idx = np.argmin(np.abs(T_vals - T_test))
    print(f"  T={T_vals[idx]:.3f} keV -> Qy={Qy_vals[idx]:.4f} e-/keV")
