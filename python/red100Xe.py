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

Ne_paper = np.array([0,1,2,3,4,5,6,7,8,9,10], dtype=float)
extraidos_paper = np.array([13.22803474907305, 7.892166913454116, 1.332326272355366,
                            0.20053061482507134, 0.03018219209250505, 0.004289421792462578,
                            0.0005756050699886731, 7.72414587856261e-05, 1.38097849568979e-05,
                            1.313348741621312e-06, 2.633666306370651e-07])
creados_paper = np.array([8.851981516921079, 14.009332579052515, 6.644004325810993,
                          2.365003434154016, 0.9442301890135518, 0.3992510858028725,
                          0.14211763609865755, 0.056740620524128794, 0.012050272964844067,
                          0.006788529364591986, 0.002416448646450225])

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
