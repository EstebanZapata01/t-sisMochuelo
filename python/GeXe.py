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
file_red = outdir + 'chi2_nsi_2DXe.dat'         # Datos de RED-100 (Xe)
file_conus = outdir + 'chi2_nsi_2Dconus.dat'  # Datos de CONUS+ (Ge)
config_file = outdir + 'nsi_configXe.txt'       # Etiquetas del caso actual

# Leer etiquetas (asumiendo que ipar es el mismo para ambos)
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

# Cargar ambos datasets
x_red, y_red, z_red = get_contour_data(file_red)
x_con, y_con, z_con = get_contour_data(file_conus)

# Nivel para 90% C.L. (2 d.o.f)
level_90 = 4.61

# ============================================================================
# GRÁFICA DE COMPARACIÓN
# ============================================================================
plt.figure(figsize=(8, 7))

# --- RED-100 (Sombreado Azul) ---
plt.contourf(x_red, y_red, z_red, levels=[0, level_90], colors=['#1f77b4'], alpha=0.3)
red_line = plt.contour(x_red, y_red, z_red, levels=[level_90], colors=['#1f77b4'], linewidths=2)

# --- CONUS+ (Línea Roja / Estilo del artículo) ---
conus_line = plt.contour(x_con, y_con, z_con, levels=[level_90], colors=['red'], 
                         linestyles='--', linewidths=2)

# Referencias
plt.axhline(0, color='black', lw=0.5, alpha=0.5)
plt.axvline(0, color='black', lw=0.5, alpha=0.5)

# Límites (ajusta según el caso, p.ej. ipar=5)
plt.xlim(-0.4, 0.8)
plt.ylim(-0.4, 0.6)

plt.xlabel(xlbl, fontsize=14)
plt.ylabel(ylbl, fontsize=14)
plt.title(f'Comparación de Sensibilidad NSI (90% C.L.)\nLXe (RED-100) vs Ge (CONUS+)', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.4)

# Leyenda manual
legend_elements = [
    Line2D([0], [0], color='#1f77b4', lw=2, label='RED-100 (Xe)'),
    Line2D([0], [0], color='red', lw=2, linestyle='--', label='CONUS+ (Ge)')
]
plt.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig(outdir + 'comparacion_nsi_xe_ge.png', dpi=300)
plt.show()
