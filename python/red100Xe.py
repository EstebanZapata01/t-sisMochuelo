#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para graficar los resultados de red100_nest.f90.
Genera:
  - fig_ionizationXe.png: Espectro en electrones de ionización (incluye datos del paper)

(El panel de retroceso Kopeikin/Mueller/combinado que este script generaba
por separado, fig_recoilXe.png, quedó cubierto por
python/recoil_spectrum_XeAr.py -> fig_recoil_XeAr.png, que compara Xe y Ar
con el mismo motor; se quitó de aquí para no duplicar salidas.)
"""

import numpy as np
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()
import os

# ===================== CONFIGURACIÓN =====================
datadir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'

# Fig. 3 del paper digitalizada (datos/red100_fig3_digitalizado.csv)
_d3 = np.loadtxt(datadir + 'red100_fig3_digitalizado.csv', delimiter=',', skiprows=1)
Ne_paper, creados_paper, extraidos_paper = _d3[:, 0], _d3[:, 1], _d3[:, 2]

# ===================== FIGURA: ELECTRONES DE IONIZACIÓN =====================
print("Generando figura: Espectro en electrones de ionización...")
file_ion = datadir + 'ionization_electrones.dat'

if os.path.exists(file_ion):
    data_i = np.loadtxt(file_ion)
    Ne = data_i[:, 0]
    creados = data_i[:, 1]
    extraidos = data_i[:, 2]
    
    fig2, ax2 = plt.subplots(figsize=(10, 6))

    # Tus datos (simulación)
    ax2.scatter(Ne, creados, c='k', marker='o', label='Creados (sim.)', s=30, zorder=3)
    ax2.scatter(Ne, extraidos, c='#8f4444', marker='s', label=f'Extraídos (sim., EEE≈0.33)', s=30, zorder=3)
    ax2.plot(Ne, creados, color='k', alpha=0.2, lw=1, zorder=2)
    ax2.plot(Ne, extraidos, color='#8f4444', alpha=0.2, lw=1, zorder=2)

    # Datos del paper RED-100
    ax2.scatter(Ne_paper, creados_paper, c='#33546e', marker='^', label='Creados (paper)', s=40, zorder=4)
    ax2.scatter(Ne_paper, extraidos_paper, c='#5c7053', marker='v', label='Extraídos (paper)', s=40, zorder=4)

    ax2.set_xlabel('Número de electrones de ionización')
    ax2.set_ylabel('Eventos / (kg · día)')
    ax2.set_title('Espectro en electrones de ionización – RED-100 (Xe)', fontsize=14)
    ax2.legend()
    ax2.set_yscale('log')
    
    ax2.set_xlim(-0.5, 10.5)
    ax2.set_xticks([0,2,4,6,8,10])
    ax2.set_ylim(1e-8, 1e2)
    ax2.set_yticks([1e-8,1e-4,1])

    fig2.tight_layout()
    fig2.savefig(datadir + 'fig_ionizationXe.png', dpi=300)
    plt.close(fig2)
    print("  → fig_ionizationXe.png guardado.")
else:
    print(f"  [ERROR] No se encontró {file_ion}")
