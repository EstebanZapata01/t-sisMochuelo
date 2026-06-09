#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re

data_dir = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
files = sorted(glob.glob(os.path.join(data_dir, "dsig_dEer_Enu*.dat")))

if not files:
    print("No se encontraron archivos con el patrón 'dsig_dEer_Enu*.dat' en", data_dir)
    exit(1)

colors = plt.cm.viridis(np.linspace(0, 1, len(files)))
plt.figure(figsize=(8, 6))

for i, filename in enumerate(files):
    # Extraer energía del neutrino permitiendo espacios opcionales
    basename = os.path.basename(filename)
    match = re.search(r'dsig_dEer_Enu\s*([\d.]+)\.dat', basename)
    if not match:
        print(f"No se pudo extraer energía de {filename}, se omite.")
        continue
    Enu = float(match.group(1))
    label = f"$E_\\nu = {Enu:.1f}$ MeV"

    data = np.loadtxt(filename)
    if data.ndim == 1 or data.shape[1] < 3:
        print(f"Archivo {filename} no tiene al menos 3 columnas, se omite.")
        continue

    Eer = data[:, 0]
    T= data[:,1]
    dsigma_lind = data[:, 2]   # con quenching
    dsigma_noQF = data[:, 3]    # sin quenching

    plt.plot(Eer, dsigma_lind, '-', color=colors[i], label=f'{label} (Lindhard)')
    plt.plot(T, dsigma_noQF, '--', color=colors[i], label=f'{label} (QF=1)')

plt.xlabel(r'$E_{\mathrm{er}}$ [MeV]', fontsize=12)
plt.ylabel(r'$\frac{d\sigma}{dE_{\mathrm{er}}}$ [cm$^2$/MeV]', fontsize=12)
plt.title('Sección eficaz diferencial de CEvNS con y sin quenching (Ge)', fontsize=13)
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=8)
plt.tight_layout()
plt.show()
