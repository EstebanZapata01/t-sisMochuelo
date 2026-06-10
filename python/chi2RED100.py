#!/usr/bin/env python3
"""
graficar_chi2.py
================
Genera dos figuras a partir de los archivos de salida de chi2_ON_OFF_1D.f90:

  Fig. 1  ->  Residuo ON-OFF con banda naranja al 90% C.L.  (Fig. 8 del articulo)
  Fig. 2  ->  Perfil Delta-chi2(A)  con linea al 90% C.L.   (Fig. 9 del articulo)

Uso:
    python3 graficar_chi2.py

Requiere:
    chi2_ON_OFF_banda.dat    (salida de chi2_ON_OFF_1D.f90)
    chi2_ON_OFF_perfil.dat   (salida de chi2_ON_OFF_1D.f90)
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

# ============================================================
# RUTAS  (ajusta solo este bloque)
# ============================================================
DATADIR  = Path('/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos')
F_BANDA  = DATADIR / 'chi2_ON_OFF_banda.dat'
F_PERFIL = DATADIR / 'chi2_ON_OFF_perfil.dat'
F_FIG8   = DATADIR / 'fig8_residuo_ON_OFF.pdf'
F_FIG9   = DATADIR / 'fig9_chi2_perfil.pdf'

# ============================================================
# FUNCION: lector robusto de .dat de Fortran
# Fortran escribe numeros muy pequenos como  0.1234-102  (sin E),
# este lector los convierte a  0.1234E-102  antes de parsear.
# ============================================================
def leer_dat(path, comentario='#'):
    meta = {}
    filas = []
    with open(path, encoding='utf-8') as f:
        for linea in f:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith(comentario):
                # Extraer metadatos del encabezado:  # clave = valor
                if '=' in linea:
                    partes = linea.lstrip('#').split('=', 1)
                    clave  = partes[0].strip()
                    try:
                        meta[clave] = float(partes[1].strip())
                    except ValueError:
                        pass
                continue
            # Corregir notacion Fortran sin E: 1.234-05 -> 1.234E-05
            linea = re.sub(r'(\d)([-+])(\d{2,3})\b', r'\1E\2\3', linea)
            try:
                filas.append([float(x) for x in linea.split()])
            except ValueError:
                pass
    return np.array(filas), meta


# ============================================================
# LEER ARCHIVOS
# ============================================================
print('Leyendo archivos de Fortran...')
banda,  meta_b = leer_dat(F_BANDA)
perfil, _      = leer_dat(F_PERFIL)

# --- banda: PE  R_SM  band_sup  band_inf  dN_data  sigma
pe       = banda[:, 0]
R_SM     = banda[:, 1]
band_up  = banda[:, 2]
# band_dn ya no se usa físicamente, ignoramos banda[:, 3]
dN       = banda[:, 4]
sigma    = banda[:, 5]

# --- perfil: A  chi2_sin  chi2_con
A_vals   = perfil[:, 0]
chi2_sin = perfil[:, 1]
chi2_con = perfil[:, 2]

# --- metadatos
A_best_sin = meta_b.get('A_best (sin nuisance)', np.nan)
A_90_sin   = meta_b.get('A_90   (sin nuisance)', np.nan)
A_best_con = meta_b.get('A_best (con nuisance)', np.nan)
A_90_con   = meta_b.get('A_90   (con nuisance)', np.nan)
chi2_min   = chi2_con.min()

print(f'  Bins de datos:              {len(pe)}')
print(f'  A_best (con nuisance) = {A_best_con:.4f}')
print(f'  A_90   (con nuisance) = {A_90_con:.4f}')
print(f'  chi2_min              = {chi2_min:.4f}')


# ============================================================
# ESTILO GLOBAL  (imita el estilo limpio del articulo)
# ============================================================
plt.rcParams.update({
    'font.family'     : 'serif',
    'font.size'       : 10,
    'axes.labelsize'  : 11,
    'axes.titlesize'  : 11,
    'xtick.labelsize' : 9,
    'ytick.labelsize' : 9,
    'legend.fontsize' : 9,
    'figure.dpi'      : 150,
    'axes.linewidth'  : 0.8,
    'xtick.direction' : 'in',
    'ytick.direction' : 'in',
    'xtick.top'       : True,
    'ytick.right'     : True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
})

NARANJA = '#E87722'   # color de la banda del articulo


# ============================================================
# FIGURA 1: Residuo ON-OFF  (Figura 8 del articulo)
# ============================================================
fig1, ax1 = plt.subplots(figsize=(6.8, 3.2))

# 1. Linea en cero (base del fondo sin señal)
ax1.axhline(0.0, color='0.55', lw=0.6, ls='--', zorder=1)

# 2. Predicción del Modelo Estándar (A=1)
ax1.step(pe, R_SM, where='mid', color='blue', lw=1.5, 
         linestyle='-', label='Predicción CEvNS (SM)', zorder=2)

# 3. Banda naranja: Limite superior al 90% C.L. (desde 0 hasta band_up)
#    Bajamos un poco el alpha para que deje ver la línea azul del SM debajo
ax1.fill_between(pe, 0, band_up,
                 step='mid',
                 color=NARANJA, alpha=0.3,
                 label='Límite 90% C.L.', zorder=3)

#    Linea de borde superior de la banda
ax1.step(pe, band_up, color=NARANJA, lw=1.5, where='mid', alpha=0.8, zorder=4)

# 4. Puntos ON-OFF con barras de error (zorder=5 los pone al frente del todo)
ax1.errorbar(pe, dN, yerr=sigma,
             fmt='ko', ms=3.8, lw=0.9,
             capsize=2.2, capthick=0.9,
             label='Delta ON-OFF',
             zorder=5)

# Decoracion
ax1.set_xlim(pe.min() - 2, pe.max() + 2)

y_lim = max(abs(dN + sigma).max(), abs(band_up).max()) * 1.45
ax1.set_ylim(-y_lim, y_lim)

ax1.set_xlabel('Corrected energy [PE]')
ax1.set_ylabel(r'Counts$\cdot$kg$^{-1}\cdot$day$^{-1}$')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(5))
ax1.yaxis.set_major_locator(ticker.MaxNLocator(5))
ax1.grid(True, ls=':', lw=0.35, color='0.75', zorder=0)

# Legend con 'ncol=2' para que quede más horizontal y limpio si quieres
leg1 = ax1.legend(framealpha=0.95, edgecolor='0.75',
                  loc='upper right', handlelength=1.4)

fig1.tight_layout()
fig1.savefig(F_FIG8, dpi=250, bbox_inches='tight')
print(f'\nFigura 8 guardada: {F_FIG8}')


# ============================================================
# FIGURA 2: Perfil Delta-chi2(A)  (Figura 9 del articulo)
# ============================================================
dchi2_sin = chi2_sin - chi2_sin.min()
dchi2_con = chi2_con - chi2_min

fig2, ax2 = plt.subplots(figsize=(5.8, 3.8))

# Perfil con nuisance (principal)
ax2.plot(A_vals, dchi2_con,
         color='k', lw=1.8,
         label=r'$\Delta\chi^2$ (con nuisance 16.9%)')

# Perfil sin nuisance (referencia, mas tenue)
ax2.plot(A_vals, dchi2_sin,
         color='0.55', lw=1.0, ls='--',
         label=r'$\Delta\chi^2$ (sin nuisance)')

# Nivel 90% C.L.
ax2.axhline(2.706, color=NARANJA, lw=1.3, ls='--',
            label=r'90% C.L.  ($\Delta\chi^2 = 2.706$)')

# Nivel 1 sigma
ax2.axhline(1.000, color='steelblue', lw=0.9, ls=':',
            label=r'$1\sigma$  ($\Delta\chi^2 = 1.000$)')

# Linea vertical: A_best y A_90
ax2.axvline(A_best_con, color='k',      lw=0.7, ls=':', alpha=0.6)
ax2.axvline(A_90_con,   color=NARANJA,  lw=0.7, ls=':', alpha=0.8)
ax2.axvline(1.0,        color='0.55',   lw=0.7, ls=':', alpha=0.5)

# Etiquetas en la parte superior del plot
ymax_plot = 9.0
ax2.annotate(f'$A_{{best}}={A_best_con:.2f}$',
             xy=(A_best_con, 0),
             xytext=(A_best_con + 0.08, ymax_plot * 0.72),
             fontsize=8.5, color='k',
             arrowprops=dict(arrowstyle='->', color='k', lw=0.7))

ax2.annotate(f'$A_{{90\\%}}={A_90_con:.2f}$',
             xy=(A_90_con, 2.706),
             xytext=(A_90_con - 2.5, ymax_plot * 0.50), # Movido a la izquierda para no tapar si el número es grande
             fontsize=8.5, color=NARANJA,
             arrowprops=dict(arrowstyle='->', color=NARANJA, lw=0.7))

# Etiqueta SM
ax2.text(1.0, ymax_plot * 0.05, 'SM\n(A=1)',
         ha='center', va='bottom', fontsize=7.5, color='0.5')

# Decoracion
ax2.set_xlim(max(-0.8, A_best_con - 1.5),
             min(A_90_con + 1.0, A_vals.max()))
ax2.set_ylim(-0.3, ymax_plot)
ax2.set_xlabel(r'Amplitud de la senal CEvNS  $A$')
ax2.set_ylabel(r'$\Delta\chi^2 = \chi^2(A) - \chi^2_{\min}$')

# Ajuste automático del espaciado para X cuando A_90 es un número mayor a 10
ax2.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=8))
ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())

ax2.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax2.grid(True, ls=':', lw=0.35, color='0.75', zorder=0)

leg2 = ax2.legend(framealpha=0.95, edgecolor='0.75',
                  loc='upper left', handlelength=1.8)

fig2.tight_layout()
fig2.savefig(F_FIG9, dpi=250, bbox_inches='tight')
print(f'Figura 9 guardada: {F_FIG9}')

# ============================================================
# RESUMEN FINAL
# ============================================================
ndof = len(pe)
print(f'\n=== RESUMEN ===')
print(f'  Bins usados:              {ndof}')
print(f'  chi2_min / ndof:          {chi2_min:.3f} / {ndof} = {chi2_min/ndof:.3f}')
print(f'  A_best (sin nuisance):    {A_best_sin:.4f}')
print(f'  A_90   (sin nuisance):    {A_90_sin:.4f}  x SM')
print(f'  A_best (con nuisance):    {A_best_con:.4f}')
print(f'  A_90   (con nuisance):    {A_90_con:.4f}  x SM')
if not np.isnan(A_best_con):
    if 0.0 <= A_best_con <= 1.5:
        print('  Interpretacion: compatible con la prediccion SM')
    elif A_best_con < 0.0 or A_best_con == 0.0:
        print('  Interpretacion: senal no requerida, solo limite superior')

plt.show()
