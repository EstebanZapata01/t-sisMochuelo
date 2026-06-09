#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import os


data_dir = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"


isotopes = ["U235", "U238", "Pu239", "Pu241"]
labels   = [r"$^{235}$U", r"$^{238}$U", r"$^{239}$Pu", r"$^{241}$Pu"]
colors   = ["green", "orange", "red", "m"]

plt.figure(figsize=(8,6))

for iso, lab, col in zip(isotopes, labels, colors):
    filepath = os.path.join(data_dir, f"spectrum_{iso}.dat")
    E, spec = np.loadtxt(filepath, unpack=True)
    plt.semilogy(E, spec, label=lab, color=col, linewidth=1.5)

plt.xlabel(r"$E_\nu$ (MeV)", fontsize=12)
plt.ylabel(r"$\lambda(E_\nu)$ [MeV$^{-1}$ fission$^{-1}$]", fontsize=12)
plt.title("Espectros de antineutrinos (Mueller et al. 2011)", fontsize=13)
plt.legend()
plt.grid(True, which="both", ls="--", lw=0.5)

plt.yscale("log")
plt.ylim(1e-5, 2e0)
plt.xlim(2, 8)

plt.tight_layout()
plt.show()

