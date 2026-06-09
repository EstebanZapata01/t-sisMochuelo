#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker

# ============================================================================
# PARÁMETROS
# ============================================================================
outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
path_data = outdir + 'chi2_sin2theta.dat'
path_plot = outdir + 'grafica_chi2.png'
path_bestfit = outdir + 'best_fit.txt'

figsize = (8, 7)
xlim = (0.01, 0.50)
ylim = (0, 10)

n_interp = 1000
kind_interp = 'cubic'

s2w_sm = 0.23857

# ============================================================================
# CARGA DE DATOS
# ============================================================================
data = np.loadtxt(path_data)
s2w = data[:, 0]
chi2 = data[:, 1]

# ============================================================================
# ANÁLISIS
# ============================================================================
idx_min = np.argmin(chi2)
s2w_min = s2w[idx_min]
chi2_min = chi2[idx_min]
delta_chi2 = chi2 - chi2_min

# Interpolación
f_interp = interp1d(s2w, delta_chi2, kind=kind_interp)
s2w_fine = np.linspace(s2w.min(), s2w.max(), n_interp)
delta_fine = f_interp(s2w_fine)

# Intervalo 1σ
indices_1s = np.where(delta_fine <= 1.0)[0]
if len(indices_1s) > 0:
    lim_inf_1s = s2w_fine[indices_1s[0]]
    lim_sup_1s = s2w_fine[indices_1s[-1]]
    err = 0.5 * ((lim_sup_1s - s2w_min) + (s2w_min - lim_inf_1s))
else:
    lim_inf_1s = lim_sup_1s = err = np.nan

# ============================================================================
# GRÁFICA
# ============================================================================
fig, ax = plt.subplots(figsize=figsize)

# Curva principal
line_data, = ax.plot(s2w, delta_chi2, color='black', lw=1.6,
                     label=r'$\Delta\chi^2$')

# Niveles de confianza
levels = [
    (1.00, r'1$\sigma$'),
    (2.71, r'90% CL'),
    (3.84, r'2$\sigma$'),
    (6.63, r'99% CL'),
    (9.00, r'3$\sigma$'),
]

for y, _ in levels:
    ax.axhline(y, color='black', lw=0.8)

# SM y Best fit
line_sm = ax.axvline(s2w_sm, color='black', linestyle=':', lw=1.2,
                     label=fr'SM = {s2w_sm:.4f}')

line_bf = ax.axvline(s2w_min, color='black', linestyle='--', lw=1.2,
                     label=fr'Best fit = {s2w_min:.4f}')

# ============================================================================
# ETIQUETAS CL FUERA DEL GRÁFICO
# ============================================================================
for y, lab in levels:
    ax.text(1.02, y, lab,
            transform=ax.get_yaxis_transform(),
            ha='left', va='center', fontsize=9)

# ============================================================================
# LEYENDA ÚNICA (TODO AQUÍ)
# ============================================================================
results_text = (
    fr'1$\sigma$: [{lim_inf_1s:.4f}, {lim_sup_1s:.4f}]' '\n'
    fr'$\pm {err:.4f}$'
)

handles = [line_data, line_sm, line_bf]
labels = [
    r'$\Delta\chi^2$',
    fr'SM = {s2w_sm:.4f}',
    fr'Best fit = {s2w_min:.4f}' + '\n' + results_text
]

ax.legend(handles, labels,
          loc='upper left',
          fontsize=9,
          frameon=True,
          edgecolor='black')

# ============================================================================
# ESTILO
# ============================================================================
ax.set_xlabel(r'$\sin^2\theta_W$', fontsize=12)
ax.set_ylabel(r'$\Delta\chi^2$', fontsize=12)

ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.xaxis.set_major_locator(ticker.MultipleLocator(0.10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.02))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax.tick_params(direction='in', which='both',
               top=True, right=True,
               length=6)

# Marco
for spine in ax.spines.values():
    spine.set_linewidth(1.0)

# Margen derecho extra para texto CL
plt.subplots_adjust(right=0.82)

plt.tight_layout()
plt.savefig(path_plot, dpi=300)
plt.show()

# ============================================================================
# SALIDA
# ============================================================================
print("\n--- RESULTADOS ---")
print(f"Best fit: {s2w_min:.6f}")
print(f"1σ: [{lim_inf_1s:.6f}, {lim_sup_1s:.6f}]")
print(f"Error: ±{err:.6f}")
print(f"SM: {s2w_sm:.6f}")

with open(path_bestfit, 'w') as f:
    f.write(f"sin^2theta_W = {s2w_min:.8f} ± {err:.8f}\n")
    f.write(f"1σ interval = [{lim_inf_1s:.8f}, {lim_sup_1s:.8f}]\n")
    f.write(f"SM = {s2w_sm:.8f}\n")
