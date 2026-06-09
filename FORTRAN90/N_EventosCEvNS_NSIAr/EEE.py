#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg') # Modo silencioso para servidores/terminales

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import AutoMinorLocator
import sys

print("==================================================")
print("  Visualizador de Sensibilidad NSI - RED-100 (Ar) ")
print("==================================================")

print("\n[1/5] Configurando entorno para alta resolución...")

# 1. Cargar los datos (Archivos de Argón)
print("[2/5] Cargando archivos .dat de Argón (1,000,000 de puntos por archivo)...")
try:
    data_pes = np.loadtxt('chi2_ar_pesimista.dat')
    data_nom = np.loadtxt('chi2_ar_nominal.dat')
    data_opt = np.loadtxt('chi2_ar_optimista.dat')
    print("      -> ¡Datos cargados con éxito!")
except OSError:
    print("\n[ERROR] No se encontraron los archivos .dat de Argón.")
    print("Asegúrate de que chi2_ar_pesimista.dat, chi2_ar_nominal.dat y chi2_ar_optimista.dat están en el directorio actual.")
    sys.exit(1)

# 2. Procesamiento de la malla (Ajustado a 1000x1000)
print("[3/5] Redimensionando matrices a 1000x1000...")
resolucion = 1000 

u_vals = data_nom[:, 0].reshape(resolucion, resolucion)
d_vals = data_nom[:, 1].reshape(resolucion, resolucion)

chi2_pes = data_pes[:, 2].reshape(resolucion, resolucion)
chi2_nom = data_nom[:, 2].reshape(resolucion, resolucion)
chi2_opt = data_opt[:, 2].reshape(resolucion, resolucion)

# 3. Configuración de la gráfica
print("[4/5] Generando contornos al 90% C.L....")
fig, ax = plt.subplots(figsize=(9, 7))

# Título adaptado para Argón
ax.set_title(r'Sensibilidad a NSI en LAr: Impacto de la EEE', fontsize=16, pad=15, fontweight='bold')
ax.set_xlabel(r'$\epsilon_{ee}^{dV}$', fontsize=14)
ax.set_ylabel(r'$\epsilon_{ee}^{uV}$', fontsize=14)

nivel_cl = [4.61] # 90% C.L.

# Dibujar las bandas de sensibilidad
ax.contour(u_vals, d_vals, chi2_pes, levels=nivel_cl, colors='#d62728', linestyles='--', linewidths=2.5)
ax.contour(u_vals, d_vals, chi2_nom, levels=nivel_cl, colors='#1f77b4', linestyles='-', linewidths=2.5)
ax.contour(u_vals, d_vals, chi2_opt, levels=nivel_cl, colors='#2ca02c', linestyles='-.', linewidths=2.5)

# Detalles estéticos
ax.axhline(0, color='black', lw=1, ls=':')
ax.axvline(0, color='black', lw=1, ls=':')
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)

# Leyenda con los valores de EEE para Argón
leg_pes = mpatches.Patch(color='#d62728', label='EEE = 0.40 (Pesimista)')
leg_nom = mpatches.Patch(color='#1f77b4', label='EEE = 0.70 (Nominal)')
leg_opt = mpatches.Patch(color='#2ca02c', label='EEE = 1.00 (Optimista)')
ax.legend(handles=[leg_opt, leg_nom, leg_pes], loc='upper right', fontsize=12, shadow=True)

plt.tight_layout()

# 5. Guardar con nombre para Argón
out_name = 'contornos_nsi_ar_eee_HD.png'
print(f"[5/5] Exportando imagen de alta resolución: {out_name}...")
plt.savefig(out_name, dpi=300, bbox_inches='tight')

print("\n==================================================")
print(f"  ¡Listo! Revisa el archivo {out_name} ")
print("==================================================")
