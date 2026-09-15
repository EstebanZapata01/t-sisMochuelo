#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sensibilidad de Xe vs. exposicion (rama esperada, Asimov de conteo puro).

RED-100 ajusta solo la amplitud: no hay nuisance de flujo (el 16.9% que usa
CONUS+ es un presupuesto sistematico especifico de germanio). Por tanto la
proyeccion estadistica

    A_90_esperado(exposicion) = 1 + sqrt(2.706 / (mult * S2))

decrece monotonamente y TIENDE A 1 cuando la exposicion -> infinito: la
estadistica sola NO impone un piso. Lo que realmente acota el alcance de
RED-100 son los sistematicos discretos (modelo de espectro 63-94, yield NEST
27-135, EEE 43-78 xSM; ver metodologia.tex, "brecha ideal->real") y la
extrapolacion a 1 anio del propio SVII del paper (15-20 xSM, escalado
efectivo ~t^-0.2). Esta figura contrasta esas dos lecturas.

Solo Xe: es el unico blanco con datos reales, asi que comparar "proyeccion
esperada" vs. "dato real 2024" vs. "cita del propio paper (SVII)" es una
comparacion dentro del MISMO experimento, siempre legitima. (El analisis de
Ar va en python/ar_fondo_validacion.py, validado solo contra el S/sqrt(B)
que declara la ref.[46], sin barrer exposicion ni compararse con Xe.)

Entrada: datos/sensib_real_Xe.dat (chi2.f90, Sec. 5d):
    # mult  exposicion_kgd  A_90_esperado

Salida: datos/fig_sensibilidad_exposicion_completa.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, GRIS, NEGRO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
EXPO_BASE = 192.0          # kg*dia, exposicion base real (volumen fiducial 2024)
A90_REAL_2024 = 107.32     # chi2.f90, ajuste ON-OFF real (sin nuisance)


def load_real_xe():
    """sensib_real_Xe.dat -> expo[kg dia], A90_esperado."""
    rows = []
    with open(f"{BASE}/sensib_real_Xe.dat", encoding="utf-8") as f:
        for ln in f:
            if ln.lstrip().startswith("#") or not ln.strip():
                continue
            p = ln.split()
            rows.append([float(p[1]), float(p[2])])
    a = np.array(rows)
    return a[:, 0], a[:, 1]

expo_xe, a90_esp = load_real_xe()

print("=" * 78)
print(" SENSIBILIDAD de Xe vs EXPOSICION -- rama esperada (Asimov, sin nuisance)")
print("=" * 78)
print(f" Xe esperado : {a90_esp[0]:.2f} (@{expo_xe[0]:.0f} kgd) -> "
      f"{a90_esp[-1]:.2f} (@{expo_xe[-1]:.0f} kgd), decreciente hacia 1")

# --------------------------------------------------------------- puntos citados
PAPER_1YR_EXPO = 200.0 * EXPO_BASE   # SVII: 1 mes OFF + 11 meses ON
PAPER_1YR_A90 = (15.0, 20.0)         # SVII, cifra publicada, no recalculada

# rangos sistematicos discretos del propio paper (base SM2018), como banda
# de referencia -- NO se escalan con la exposicion
SYS_SPECTRO = (63.0, 94.0)   # 4 modelos de espectro, Tabla I
SYS_NEST    = (27.0, 135.0)  # yield NEST a los bordes de banda, SVII
SYS_EEE     = (43.0, 78.0)   # EEE 32.8 +/- 2.8%, SVII

# --------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=(7.2, 5.0))

# banda de sistematicos discretos (rango envolvente), franja horizontal
sys_lo = min(SYS_SPECTRO[0], SYS_NEST[0], SYS_EEE[0])
sys_hi = max(SYS_SPECTRO[1], SYS_NEST[1], SYS_EEE[1])
ax.axhspan(sys_lo, sys_hi, color=GRIS, alpha=0.16, lw=0, zorder=0,
           label=f"sistemáticos discretos ({sys_lo:.0f}–{sys_hi:.0f}×SM)")

ax.plot(expo_xe, a90_esp, "o-", color=C_XE, ms=4,
        label="Xe, proyección esperada (Asimov)")

ax.plot([EXPO_BASE], [A90_REAL_2024], "D", color=NEGRO, ms=7, mfc="white", mew=1.3,
        label="Xe, dato real 2024")
ax.errorbar([PAPER_1YR_EXPO], [np.mean(PAPER_1YR_A90)],
            yerr=[[np.mean(PAPER_1YR_A90) - PAPER_1YR_A90[0]],
                  [PAPER_1YR_A90[1] - np.mean(PAPER_1YR_A90)]],
            fmt="^", color=GRIS, ms=8, capsize=3, lw=1.1,
            label="extrapolación a 1 año, §VII")

ax.axhline(1.0, color=C_XE, lw=0.9, ls=":")
ax.text(expo_xe[1], 1.03, r"límite estadístico ($A_{90}\to 1$)",
        fontsize=8.0, va="bottom", color=C_XE)

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Exposición reactor ON [kg·día]")
ax.set_ylabel(r"$A_{90}$ [$\times$SM]")
ax.legend(loc="upper right")
ax.set_ylim(0.9, 200)

fig.tight_layout()
fig.savefig(f"{BASE}/fig_sensibilidad_exposicion_completa.png")
print(f"\n  {BASE}/fig_sensibilidad_exposicion_completa.png")
