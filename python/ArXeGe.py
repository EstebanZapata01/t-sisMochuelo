#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

# ============================================================================
# CONFIGURACIÓN DE RUTAS
# ============================================================================
outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
file_xe = outdir + 'chi2_nsi_2DXe.dat'         # Datos de RED-100 (LXe)
file_ar = outdir + 'chi2_nsi_2DAr.dat'         # Datos de RED-100 (LAr)
file_ge = outdir + 'chi2_nsi_2Dconus.dat'      # Datos de CONUS+ (Ge)
config_file = outdir + 'nsi_configAr.txt'      # Usamos el config más reciente

# Leer etiquetas
if os.path.exists(config_file):
    with open(config_file, 'r') as f:
        lines = f.readlines()
        xlbl = lines[0].strip()
        ylbl = lines[1].strip()
else:
    xlbl, ylbl = r'$\epsilon_{x}$', r'$\epsilon_{y}$'

# ============================================================================
# CARGA Y PROCESAMIENTO
# ============================================================================
def get_contour_data(filename):
    data = np.loadtxt(filename)
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    # Reshape (ny, nx) asumiendo malla de 1000x1000
    z = data[:, 2].reshape(len(y), len(x))
    z = z - np.min(z) # Delta Chi2
    return x, y, z

# Cargar los tres datasets
x_xe, y_xe, z_xe = get_contour_data(file_xe)
x_ar, y_ar, z_ar = get_contour_data(file_ar)
x_ge, y_ge, z_ge = get_contour_data(file_ge)

# Nivel para 90% C.L. (2 d.o.f)
level_90 = 4.61

# ============================================================================
# GRÁFICA DE COMPARACIÓN MÚLTIPLE
# ============================================================================
plt.figure(figsize=(8, 7))

# --- RED-100 LXe (Azul - Base) ---
plt.contourf(x_xe, y_xe, z_xe, levels=[0, level_90], colors=['#1f77b4'], alpha=0.3)
plt.contour(x_xe, y_xe, z_xe, levels=[level_90], colors=['#1f77b4'], linewidths=2)

# --- RED-100 LAr (Verde - Intersección) ---
plt.contourf(x_ar, y_ar, z_ar, levels=[0, level_90], colors=['#2ca02c'], alpha=0.15)
plt.contour(x_ar, y_ar, z_ar, levels=[level_90], colors=['#2ca02c'], linewidths=2)

# --- CONUS+ Ge (Rojo - Punteado, estilo literatura) ---
plt.contour(x_ge, y_ge, z_ge, levels=[level_90], colors=['red'], linestyles='--', linewidths=2)

# Referencias cruzadas en cero
plt.axhline(0, color='black', lw=0.5, alpha=0.5)
plt.axvline(0, color='black', lw=0.5, alpha=0.5)

# Límites (ajustables según la vista que necesites)
plt.xlim(-0.4, 0.8)
plt.ylim(-0.4, 0.6)

plt.xlabel(xlbl, fontsize=14)
plt.ylabel(ylbl, fontsize=14)
plt.title('Comparación de Sensibilidad NSI (90% C.L.)\nLXe vs LAr vs Ge', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.4)

# Leyenda manual combinada
legend_elements = [
    Line2D([0], [0], color='#1f77b4', lw=2, label='RED-100 (LXe)'),
    Line2D([0], [0], color='#2ca02c', lw=2, label='RED-100 (LAr)'),
    Line2D([0], [0], color='red', lw=2, linestyle='--', label='CONUS+ (Ge)')
]
plt.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig(outdir + 'comparacion_nsi_Xe_Ar_Ge.png', dpi=300)
plt.show()
