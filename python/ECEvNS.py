#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os

# Ruta del archivo generado por Fortran (eventos por kg)
filepath = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/eventos_conus_total.dat'
data = np.loadtxt(filepath)
bin_centers = data[:, 0]      # centros de los bins en eV
eventos_kg = data[:, 1]       # predicción SM en kg⁻¹ (10 eV)⁻¹

# Datos experimentales de CONUS+ (directamente del gráfico, en las mismas unidades)
exp_counts = np.array([13.27, 38.06, 29.09, 21.37, 12.27, 4.92, 4.29, 0.06, 9.15,
                       2.42, 3.92, -2.55, 0.56, -10.52, -4.29, -8.78, -0.56, -0.93, 13.52])
exp_errors = np.array([22.99, 15.70, 12.64, 12.21, 11.52, 11.58, 11.27, 10.77, 11.02,
                       10.28, 10.34, 10.77, 10.34, 10.15, 10.09, 9.90, 9.84, 9.78, 9.59])

# Graficar
plt.figure(figsize=(10,6))
plt.bar(bin_centers, eventos_kg, width=9, align='center', alpha=0.7, color='blue',
        label='Predicción SM (este trabajo)')
plt.errorbar(bin_centers, exp_counts, yerr=exp_errors, fmt='o', color='black',
             capsize=3, label='Datos CONUS+ (excesos)')
plt.xlabel('Energía reconstruida (eV$_{ee}$)', fontsize=12)
plt.ylabel('Excess counts [kg$^{-1}$ (10 eV)$^{-1}$]', fontsize=12)
plt.title('CEvNS en CONUS+', fontsize=14)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(150, 360)
plt.tight_layout()
plt.savefig('comparacion_CONUS+.png', dpi=300)
print("Gráfico guardado como 'comparacion_CONUS+.png'")
plt.show()
