#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sensibilidad de Xe vs. exposicion ON, en la CONFIGURACION REAL (residuo
ON-OFF, chi2.f90), no en el limite estadistico idealizado.

El punto es que la resta ON-OFF tiene DOS presupuestos de Poisson (ON y
OFF), y la extrapolacion publicada a un ano (SVII: 1 mes OFF + 11 meses ON)
da solo 15-20xSM porque la muestra OFF (que fija el fondo) NO crece igual de
rapido que la ON. Se comparan tres supuestos sobre como crece t_OFF cuando
t_ON crece:
  (a) t_OFF/t_ON fijo = 61/192 (la proporcion fiduciaria de 2024: toda la
      varianza escala con la exposicion, como en la version anterior de esta
      figura);
  (b) t_OFF/t_ON = 1/11 (un mes apagado por cada once meses encendido, el
      supuesto de SVII);
  (c) t_OFF fijo en 61 kg dia (el mismo dato OFF de 2024, solo crece t_ON):
      el piso que impone una muestra OFF que no vuelve a tomarse.

Las tasas de fondo R_ON,k y R_OFF,k por bin de energia corregida [PE] se
RECONSTRUYEN del residuo ON-OFF publicado y digitalizado (delta_ON_OFF,
sigma_stat en datos/chi2_ON_OFF_banda.dat) resolviendo, para cada bin,
    R_ON - R_OFF = delta_ON_OFF ,   R_ON/t_ON0 + R_OFF/t_OFF0 = sigma_stat^2 ,
con t_ON0=192, t_OFF0=61 kg dia (exposiciones FIDUCIARIAS de 2024,
Codigos/doc/papers_resumen_es.tex l.686-687), suponiendo sigma_stat
puramente de Poisson (las barras de error del residuo publicado). Los bines
que salen negativos por ruido de digitalizacion se capan a 0 (se cuentan en
stdout). Con R_ON,k, R_OFF,k fijos (tasas, no cuentas) se recalcula
sigma_k^2(t_ON) = R_ON,k/t_ON + R_OFF,k/t_OFF(t_ON) para cada escenario y
    A_90(t_ON) = 1 + sqrt(2.706 / Sum_k R_SM,k^2/sigma_k^2(t_ON)).

Entrada: datos/chi2_ON_OFF_banda.dat (PE_center, R_SM, delta_ON_OFF, sigma_stat)
Salida : datos/fig_sensibilidad_exposicion_completa.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, GRIS, NEGRO, FIG15
aplicar(grande=True)

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
T_ON0, T_OFF0 = 192.0, 61.0          # kg dia, exposiciones FIDUCIARIAS 2024
DCHI2_90 = 2.706

d = np.loadtxt(f"{BASE}/chi2_ON_OFF_banda.dat", comments="#")
RSM, dN, sig = d[:, 1], d[:, 4], d[:, 5]

R_OFF = (sig**2 - dN / T_ON0) / (1.0 / T_ON0 + 1.0 / T_OFF0)
R_ON = dN + R_OFF
n_neg = int((R_OFF < 0).sum() + (R_ON < 0).sum())
R_OFF = np.maximum(R_OFF, 0.0); R_ON = np.maximum(R_ON, 0.0)


def A90(t_on, escenario):
    if escenario == "a":
        t_off = t_on * (T_OFF0 / T_ON0)
    elif escenario == "b":
        t_off = t_on / 11.0
    else:
        t_off = T_OFF0
    sig2 = R_ON / t_on + R_OFF / max(t_off, 1e-9)
    S2 = np.sum(RSM**2 / sig2)
    return 1.0 + np.sqrt(DCHI2_90 / S2)


expo = np.logspace(np.log10(T_ON0), np.log10(2.0e5), 200)
CURVAS = {"a": ("proporción 2024 fija (192:61)", "-", C_XE),
          "b": ("11 ON : 1 OFF (un año de operación)", "--", "#7a5c9e"),
          "c": ("OFF fijo (61 kg·día)", "-.", "#a86a43")}

# Exposicion ON a 1 anio: 11 meses encendido de 12 (SVII), con el mismo
# cociente livetime/tiempo-real que 2024 (arXiv:2411.18641 Sec. II.B,
# "average livetime over elapsed real time ratio with CEvNS trigger is
# about 60%") sobre los T_VIVO_ON0=2.63 dias vivos que dan T_ON0=192 kg dia.
# F_VIVO es un supuesto (el articulo no da la exposicion real de ese anio).
T_VIVO_ON0 = 2.63    # dias vivos de 2024 (arXiv:2411.18641 Sec. II.B)
F_VIVO = 0.60        # fraccion livetime/tiempo-real, igual que 2024
DIAS_ON_1YR = 365.25 * 11.0 / 12.0   # 11 de 12 meses con el reactor ON (SVII)
EXPO_1YR = T_ON0 * (DIAS_ON_1YR * F_VIVO) / T_VIVO_ON0
PAPER_1YR_A90 = (15.0, 20.0)
SYS_NEST_2024 = (27.0, 135.0)   # banda de Q_y en el limite de 2024 (rendimiento NEST)
A90_OBS_2024 = 76.6204          # chi2_ON_OFF_banda.dat, cabecera "A_90" (ajuste observado real)

fig, ax = plt.subplots(figsize=FIG15)
ax.plot([T_ON0, T_ON0], SYS_NEST_2024, color=GRIS, lw=3.2, alpha=0.55, zorder=0, solid_capstyle="butt",
        label=rf"límite de 2024, $Q_y$ en su banda ({SYS_NEST_2024[0]:.0f}–{SYS_NEST_2024[1]:.0f}×SM)")

vals_1yr = {}
for k, (lab, ls, col) in CURVAS.items():
    curva = np.array([A90(t, k) for t in expo])
    ax.plot(expo, curva, ls, color=col, lw=2.0, label=lab)
    vals_1yr[k] = A90(EXPO_1YR, k)

ax.plot([T_ON0], [A90_OBS_2024], "D", color=NEGRO, ms=8, mfc="white", mew=1.4,
        label="dato real 2024 (observado)")
ax.errorbar([EXPO_1YR], [np.mean(PAPER_1YR_A90)],
            yerr=[[np.mean(PAPER_1YR_A90) - PAPER_1YR_A90[0]],
                  [PAPER_1YR_A90[1] - np.mean(PAPER_1YR_A90)]],
            fmt="^", color=GRIS, ms=9, capsize=3, lw=1.2,
            label="extrapolación a un año (publicada)")

ax.axhline(1.0, color="0.7", lw=0.9, ls=":")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Exposición reactor ON [kg·día]")
ax.set_ylabel(r"$A_{90}$ [$\times$SM]")
ax.legend(loc="lower left", fontsize=9.6)
ax.set_ylim(0.9, 200)

fig.tight_layout()
out = f"{BASE}/fig_sensibilidad_exposicion_completa.png"
fig.savefig(out, dpi=300, bbox_inches="tight")

print(f"  bines capados a 0 por reconstruccion negativa: {n_neg} / {2*len(RSM)}")
print(f"  A_90 a t_ON=192 (control, debe dar 73.07=esperado Asimov): a={A90(T_ON0,'a'):.4f} "
      f"b={A90(T_ON0,'b'):.4f} c={A90(T_ON0,'c'):.4f}")
print(f"  Exposicion 1 anio = {EXPO_1YR:.0f} kg dia "
      f"({EXPO_1YR/T_ON0:.1f}x base; SVII 11/12 meses ON, F_VIVO={F_VIVO:.2f})")
print(f"  A_90 a 1 anio: (a) proporcional={vals_1yr['a']:.3f}  (b) 11:1={vals_1yr['b']:.3f}  "
      f"(c) OFF fijo={vals_1yr['c']:.3f}")
print(f"  (c) asintotico (t_ON->infinito): {A90(1e9,'c'):.3f}")
print(f"  -> {out}")
