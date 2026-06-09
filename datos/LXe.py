#!/usr/bin/env python3
"""
Script v3: Búsqueda dinámica del mejor modelo de detector en NEST
para simular RED-100 a 218 V/cm.
"""
import nestpy
import numpy as np

print(f"nestpy versión: {nestpy.__nest_version__}")
print("Escaneando detectores disponibles y evaluando Q_y(1 keV)...\n")

# Extraer automáticamente todos los detectores compilados en tu versión de nestpy
detector_names = [d for d in dir(nestpy.detectors) if "Detector" in d and not d.startswith("__")]

# Si por la versión de pybind11 no expone el dir() correctamente, usamos una lista segura:
if not detector_names:
    detector_names = [
        "DetectorExample_LUX", "DetectorExample_XENON10", 
        "DetectorExample_XENON1T", "DetectorExample_Run03", "DetectorExample_Run04"
    ]

# Parámetros objetivo (del paper de NEST / RED-100)
TARGET_ENERGY = 1.0  # keV
TARGET_QY = 6.2      # e-/keV (Valor esperado aproximado a 218 V/cm)
TOLERANCIA = 0.6     # Rango aceptable para darle el ✓

# Parámetros estándar de la interacción NR
nr_params = [11.0, 1.1, 0.048, -0.0533, 12.6, 0.3, 2.0, 0.3, 2.0, 0.5, 1.0, 1.0]
er_params = [-1.0] * 10

best_det_name = None
best_diff = float('inf')
best_qy = 0.0

# --- 1. FASE DE PRUEBA Y SELECCIÓN ---
print(f"{'Modelo de Detector':<30} | {'Q_y(1 keV) [e-/keV]':<20} | {'Status'}")
print("-" * 65)

for name in detector_names:
    try:
        # Instanciar detector dinámicamente
        det_class = getattr(nestpy.detectors, name)
        det = det_class()
        nc = nestpy.NESTcalc(det)
        
        # Calcular yield a 1 keV
        y = nc.GetYields(interaction = nestpy.interactions.NR,
                         energy = TARGET_ENERGY,
                         density = 2.9,
                         drift_field = 218.0,
                         A = 131,
                         Z = 54,
                         nr_parameters = nr_params,
                         er_parameters = er_params)
        
        qy = y.ElectronYield / TARGET_ENERGY
        diff = abs(qy - TARGET_QY)
        
        # Evaluar si pasa el control de calidad
        status = "✓ Aprobado" if diff <= TOLERANCIA else "✗ Desviado"
        print(f"{name:<30} | {qy:<20.4f} | {status}")
        
        # Guardar el que más se acerque al objetivo
        if diff < best_diff:
            best_diff = diff
            best_det_name = name
            best_qy = qy
            
    except Exception as e:
        # Algunos detectores pueden no soportar NR o requerir otros inits
        pass

if not best_det_name:
    print("Error: No se encontró ningún detector compatible.")
    exit()

print("-" * 65)
print(f"🥇 Detector ganador: {best_det_name} (Desviación: {best_diff:.3f})\n")


# --- 2. FASE DE GENERACIÓN DE LA MATRIZ DENSA ---
print(f"Generando 'nest_218V_dense.txt' usando {best_det_name}...")

# Instanciar el mejor detector para la producción final
final_det_class = getattr(nestpy.detectors, best_det_name)
final_nc = nestpy.NESTcalc(final_det_class())

T_vals = np.linspace(0.05, 5.0, 2000)
Qy_vals = np.zeros(len(T_vals))

for i, T in enumerate(T_vals):
    try:
        y = final_nc.GetYields(interaction = nestpy.interactions.NR,
                               energy = T,
                               density = 2.9,
                               drift_field = 218.0,
                               A = 131,
                               Z = 54,
                               nr_parameters = nr_params,
                               er_parameters = er_params)
        Qy_vals[i] = y.ElectronYield / T
    except:
        Qy_vals[i] = 0.0

# Guardar en el formato que espera tu Fortran/Python
output = "nest_218V_dense.txt"
with open(output, "w") as f:
    f.write("# Q_y tabla NEST para RED-100\n")
    f.write(f"# Generado dinámicamente usando: {best_det_name}\n")
    f.write("# Xe NR, rho=2.9 g/cm3, E=218 V/cm\n")
    f.write("# T_nr[keV]   Qy[e-/keV]\n")
    f.write(f"# N_puntos = {len(T_vals)}\n")
    for T, Qy in zip(T_vals, Qy_vals):
        f.write(f"{T:.6f}   {Qy:.6f}\n")

print(f"✅ Matriz densa guardada con éxito.")
