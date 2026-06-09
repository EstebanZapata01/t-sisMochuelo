#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import matplotlib.ticker as ticker
data_dir = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
files = sorted(glob.glob(os.path.join(data_dir, "spectrum_E*.dat")))
colors = plt.cm.viridis(np.linspace(0, 1, len(files)))
plt.figure(figsize=(8,6))
for i, filename in enumerate(files):
    T, dsigma_FF1, dsigma_Helm = np.loadtxt(filename, unpack=True)
    energy_label = os.path.basename(filename).replace("spectrum_E","").replace(".dat","") + " MeV"
    plt.plot(T, dsigma_FF1, "--", color=colors[i], label=f"Eν={energy_label}, FF=1")
    plt.plot(T, dsigma_Helm, "-", color=colors[i], label=f"Eν={energy_label}, FF=Helm")
plt.xlabel("Energía de retroceso T (MeV)")
plt.ylabel("dσ/dT (cm² / MeV)")
plt.title("Sección eficaz diferencial CEvNS con y sin FF de Helm")
plt.legend(fontsize=8)
ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=15))
ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=15))
plt.grid(True, linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()
