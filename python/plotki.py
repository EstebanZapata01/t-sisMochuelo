#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pathlib import Path

path = Path("/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/K_all.dat")
data = np.loadtxt(path)
TN = data[:, 0] 
K = data[:, 1:] 

pathQF = Path("/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/K_allQF.dat")
dataQF = np.loadtxt(pathQF)
EerQF = dataQF[:, 0]
KQF = dataQF[:, 1:]

nbins = K.shape[1]

print(f"Archivo cargado: {path.name}")
print(f"Archivo cargado: {pathQF.name}")
print(f"Número de bins: {nbins}, puntos TN: {len(TN)}")

plt.figure(figsize=(9, 6))
colors = plt.cm.viridis(np.linspace(0, 1, nbins))
colorsQF = plt.cm.hsv(np.linspace(0, 1, nbins))

for i in range(nbins):
    plt.plot(TN, K[:, i],'.',
             color=colors[i],
             label=f"Bin {i+1}")

plt.xlabel(r"$E_{er}(T_N)$ [MeV]", fontsize=12)
plt.ylabel(r"$\int G(E_{er},E_{reco})$", fontsize=12)
plt.title("Función de resolución integrada por bin", fontsize=13)

for i in range(nbins):
    plt.plot(EerQF, KQF[:, i],'.',
             color=colorsQF[i],
             label=f"BinQF {i+1}")

print("TN range:", TN.min(), TN.max())
print("EerQF range:", EerQF.min(), EerQF.max())
plt.xlim(1e-4, 8e-4)
plt.ylim(0, 0.2)
#plt.yscale("log") 
plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=15))
plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=10))
plt.grid(True, ls='--', lw=0.5, alpha=0.6)
plt.legend(fontsize=8, frameon=False)
plt.tight_layout()

out_path = path.with_name("K_all_plot.png")
plt.savefig(out_path, dpi=300)
plt.show()

print(f"Gráfico guardado en: {out_path}")
