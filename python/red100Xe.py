#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para graficar los resultados de red100_nest.f90 (Versión Definitiva)
Genera:
  - fig_recoilXe.png    : Espectro de retrocesos nucleares comparativo (escala log)
  - fig_ionizationXe.png: Espectro en electrones de ionización (incluye datos del paper)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ===================== CONFIGURACIÓN =====================
datadir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'

plt.rcParams.update({'font.size': 12, 'legend.fontsize': 11, 'font.family': 'serif'})

# ===================== DATOS DEL PAPER RED-100 =====================
# (N_e, extraídos, creados) – arXiv:2403.12645
Ne_paper = np.array([0,1,2,3,4,5,6,7,8,9,10], dtype=float)
extraidos_paper = np.array([13.22803474907305, 7.892166913454116, 1.332326272355366,
                            0.20053061482507134, 0.03018219209250505, 0.004289421792462578,
                            0.0005756050699886731, 7.72414587856261e-05, 1.38097849568979e-05,
                            1.313348741621312e-06, 2.633666306370651e-07])
creados_paper = np.array([8.851981516921079, 14.009332579052515, 6.644004325810993,
                          2.365003434154016, 0.9442301890135518, 0.3992510858028725,
                          0.14211763609865755, 0.056740620524128794, 0.012050272964844067,
                          0.006788529364591986, 0.002416448646450225])

# ===================== FIGURA 1: ESPECTRO DE RETROCESOS (LOG) =====================
print("Generando Figura 1: Espectro de retrocesos nucleares...")
file_recoil = datadir + 'espectro_continuoXe.dat'

if os.path.exists(file_recoil):
    data_r = np.loadtxt(file_recoil)
    Tnr = data_r[:, 0]          # keV
    Kop = data_r[:, 1]          # Tasa Kopeikin
    Mue = data_r[:, 2]          # Tasa Mueller
    Comb = data_r[:, 3]         # Tasa Combinada

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(Tnr, Kop, label='Kopeikin completo (0-9 MeV)', lw=2)
    ax1.plot(Tnr, Mue, label='Mueller puro (2-10 MeV, KNPP)', lw=2)
    ax1.plot(Tnr, Comb, label='Combinado (Kop < 2, Mue >= 2)', lw=2, linestyle='--')
    
    ax1.set_xlabel('Energía de retroceso nuclear, $T_{nr}$ (keV)')
    ax1.set_ylabel('Eventos / (keV · kg · día)')
    ax1.set_title('Espectros de retroceso para RED-100 (Xe)', fontsize=14)
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    ax1.set_xlim(0.0, 3.5)
    ax1.set_ylim(1e-9,1.5e3)
    ax1.set_xticks([1,2,3])
    ax1.set_yticks([1e-9,1e-6,1e-3,1,1e+3])
    fig1.tight_layout()
    fig1.savefig(datadir + 'fig_recoilXe.png', dpi=300)
    plt.close(fig1)
    print("  → fig_recoilXe.png guardado.")
else:
    print(f"  [ERROR] No se encontró {file_recoil}")

# ===================== FIGURA 2: ELECTRONES DE IONIZACIÓN =====================
print("Generando Figura 2: Espectro en electrones de ionización...")
file_ion = datadir + 'ionization_electrones.dat'

if os.path.exists(file_ion):
    data_i = np.loadtxt(file_ion)
    Ne = data_i[:, 0]
    creados = data_i[:, 1]
    extraidos = data_i[:, 2]
    
    fig2, ax2 = plt.subplots(figsize=(10, 6))

    # Tus datos (simulación)
    ax2.scatter(Ne, creados, c='k', marker='o', label='Creados (sim.)', s=30, zorder=3)
    ax2.scatter(Ne, extraidos, c='r', marker='s', label=f'Extraídos (sim., EEE≈0.33)', s=30, zorder=3)
    ax2.plot(Ne, creados, color='k', alpha=0.2, lw=1, zorder=2)
    ax2.plot(Ne, extraidos, color='r', alpha=0.2, lw=1, zorder=2)

    # Datos del paper RED-100
    ax2.scatter(Ne_paper, creados_paper, c='blue', marker='^', label='Creados (paper)', s=40, zorder=4)
    ax2.scatter(Ne_paper, extraidos_paper, c='green', marker='v', label='Extraídos (paper)', s=40, zorder=4)

    ax2.set_xlabel('Número de electrones de ionización')
    ax2.set_ylabel('Eventos / (kg · día)')
    ax2.set_title('Espectro en electrones de ionización – RED-100 (Xe)', fontsize=14)
    ax2.legend()
    ax2.set_yscale('log')
    ax2.grid(True, linestyle=':', alpha=0.6)
    
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
