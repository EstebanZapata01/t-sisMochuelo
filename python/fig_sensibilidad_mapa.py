#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapa de sensibilidad Xe/Ar vs exposicion, A_90-1 y exposicion en log-log.
Reemplaza fig_sensibilidad_exposicion_completa.py, 0,9 del ancho.

Mensaje: las ideales (sin fondo) son rectas de pendiente -1/2 (A_90-1 =
sqrt(2,706/(R E))), separadas horizontalmente por el cociente de tasas
R_Ar/R_Xe (255x en N_e>=4); el fondo del Ar (B prop. a E) tambien da una
recta -1/2 pero mas arriba, y la aplana si se perfila sobre su
normalizacion (sigma_B libre, mismo S2_eff(E)=E*(S2_rate-P_rate^2/Q_rate),
tambien lineal en E).

Color = blanco (Xe azul, Ar cobre, SIEMPRE). Estilo = papel de la curva:
  punteado fino   = ideal N_e>=1
  trazos fino     = ideal N_e>=4
  continuo grueso = escenario base   (Xe: ON-OFF prop. 2024; Ar: fondo, norm. conocida)
  raya-punto grueso = variante pesimista (Xe: OFF fijo; Ar: fondo, norm. libre)

Xe real (ON-OFF, chi2.f90): ver sensibilidad_exposicion_completa.py para
la deduccion completa de R_ON, R_OFF.

Entrada: datos/chi2_ON_OFF_banda.dat, datos/sensib_ideal_{Xe,Ar}.dat,
         datos/bkg_bins_Ar.dat
Salida : datos/fig_sensibilidad_mapa.png (+ copia en figuras_metodologia)
"""
import shutil
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from estilo_tesis import aplicar, C_XE, C_AR, GRIS, NEGRO, formatear_coma, nota
aplicar()
mpl.rcParams.update({
    "axes.labelsize": 9.0, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "legend.fontsize": 7.3, "lines.linewidth": 1.15,
})

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
DEST = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Proyecto_LaTeX/figuras_metodologia"
T_ON0, T_OFF0 = 192.0, 61.0
DCHI2_90 = 2.706
LW_FINO, LW_GRUESO = 1.1, 1.7

# ------------------------------------------------- Xe real ON-OFF (chi2.f90)
banda = np.loadtxt(f"{BASE}/chi2_ON_OFF_banda.dat", comments="#")
RSM, dN, sig = banda[:, 1], banda[:, 4], banda[:, 5]
R_OFF = (sig**2 - dN / T_ON0) / (1.0 / T_ON0 + 1.0 / T_OFF0)
R_ON = dN + R_OFF
R_OFF = np.maximum(R_OFF, 0.0); R_ON = np.maximum(R_ON, 0.0)


def A90_real(t_on, t_off):
    sig2 = R_ON / t_on + R_OFF / max(t_off, 1e-9)
    S2 = np.sum(RSM**2 / sig2)
    return 1.0 + np.sqrt(DCHI2_90 / S2)


A90_OBS_2024 = 79.1394   # chi2_ON_OFF_banda.dat, cabecera "A_90"
EXPO_1YR = 14666.0   # t_ON de 1 ano = 76,4x base (cifras_tesis_helm.md #8)
PAPER_1YR = (15.0, 20.0)
A90_1YR_PROPIO = A90_real(EXPO_1YR, EXPO_1YR / 11.0) - 1.0   # OFF = ON/11

# ------------------------------------------------- ideal (sensib_ideal)
R_ideal = {}
for tag, k in (("Xe", 1), ("Xe", 4), ("Ar", 1), ("Ar", 4)):
    d = np.loadtxt(f"{BASE}/sensib_ideal_{tag}.dat", comments="#", dtype=str)
    fila = d[(d[:, 0] == str(k)) & (d[:, 1] == "Fnest")][0]
    R_ideal[(tag, k)] = float(fila[2])   # R_tot_ROI

# ------------------------------------------------- Ar con fondo (B prop. a E)
Ne, Sk, Bk0 = np.loadtxt(f"{BASE}/bkg_bins_Ar.dat", comments="#", unpack=True)
sig2b = Sk + Bk0
S2_rate_conocida = np.sum(Sk**2 / sig2b)
P_rate = np.sum(Sk * Bk0 / sig2b)
Q_rate = np.sum(Bk0**2 / sig2b)
S2_rate_libre = S2_rate_conocida - P_rate**2 / Q_rate

E = np.logspace(np.log10(10.0), np.log10(2.0e5), 300)

fig, ax = plt.subplots(figsize=(5.53, 4.8))

# ideales (finas)
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    ax.plot(E, np.sqrt(DCHI2_90 / (R_ideal[(tag, 1)] * E)), ":", color=col, lw=LW_FINO)
    ax.plot(E, np.sqrt(DCHI2_90 / (R_ideal[(tag, 4)] * E)), "--", color=col, lw=LW_FINO)

# Xe real: base (prop. 2024) y pesimista (OFF fijo)
Ereal = np.logspace(np.log10(T_ON0), np.log10(2.0e5), 200)
curva_base = np.array([A90_real(t, t * T_OFF0 / T_ON0) - 1.0 for t in Ereal])
curva_pes = np.array([A90_real(t, T_OFF0) - 1.0 for t in Ereal])
ax.plot(Ereal, curva_base, "-", color=C_XE, lw=LW_GRUESO)
ax.plot(Ereal, curva_pes, "-.", color=C_XE, lw=LW_GRUESO)

# Ar + fondo: base (conocida) y pesimista (libre)
ax.plot(E, np.sqrt(DCHI2_90 / (S2_rate_conocida * E)), "-", color=C_AR, lw=LW_GRUESO)
ax.plot(E, np.sqrt(DCHI2_90 / (S2_rate_libre * E)), "-.", color=C_AR, lw=LW_GRUESO)

# ------------------------------------------------- marcadores (sin texto en el eje)
ax.plot([T_ON0], [A90_OBS_2024 - 1.0], "D", color=NEGRO, ms=6, mfc="white", mew=1.3, zorder=6)
ax.plot([EXPO_1YR, EXPO_1YR], PAPER_1YR, "-", color="0.4", lw=7.0, alpha=0.75, solid_capstyle="butt", zorder=3)
ax.plot([EXPO_1YR], [A90_1YR_PROPIO], "o", color=NEGRO, ms=6, mfc="white", mew=1.3, zorder=6)

# ------------------------------------------------- leyenda, 3 columnas (Xe, Ar, puntos)
proxy = lambda col, ls, lw: Line2D([], [], color=col, ls=ls, lw=lw)
marker = lambda mk: Line2D([], [], color=NEGRO, marker=mk, ms=6, mfc="white", mew=1.3, ls="")
handles = [
    proxy(C_XE, ":", LW_FINO), proxy(C_XE, "--", LW_FINO),
    proxy(C_XE, "-", LW_GRUESO), proxy(C_XE, "-.", LW_GRUESO),
    proxy(C_AR, ":", LW_FINO), proxy(C_AR, "--", LW_FINO),
    proxy(C_AR, "-", LW_GRUESO), proxy(C_AR, "-.", LW_GRUESO),
    marker("D"), marker("o"), Line2D([], [], color="0.4", lw=7.0, alpha=0.75),
]
labels = [
    r"Xe ideal, $N_e\geq1$", r"Xe ideal, $N_e\geq4$",
    "Xe real, prop. 2024 (base)", "Xe real, OFF fijo (pesimista)",
    r"Ar ideal, $N_e\geq1$", r"Ar ideal, $N_e\geq4$",
    "Ar + fondo, conocida (base)", "Ar + fondo, libre (pesimista)",
    "dato real Xe 2024", "1 año, este trabajo (ON:OFF = 11:1)",
    "1 año, publicado (15–20×SM)",
]
leg = ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.17),
                 ncol=3, fontsize=6.8, frameon=False, columnspacing=1.4)

for xv, lab in ((62.0, "62"), (192.0, "192")):
    ax.axvline(xv, color="0.75", lw=0.8, ls=":")
    ax.text(xv, 0.03, lab, rotation=90, ha="right", va="bottom", fontsize=6.6, color="0.45",
            transform=ax.get_xaxis_transform())

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_ylim(4e-4, 250.0)
ax.set_xlabel("Exposición [kg·día]")
ax.set_ylabel(r"$A_{90}-1$  [$\times$SM]")
formatear_coma(ax)

fig.tight_layout()
out = f"{BASE}/fig_sensibilidad_mapa.png"
fig.savefig(out, dpi=300, bbox_extra_artists=(leg,), bbox_inches="tight")
plt.close(fig)
dest = f"{DEST}/fig_sensibilidad_mapa.png"
shutil.copyfile(out, dest)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(out) == md5(dest)

print(f"  EXPO_1YR = {EXPO_1YR:.1f} kg dia ; A90(1yr, ON:OFF=11:1) = {A90_1YR_PROPIO+1.0:.4f}"
      f" -> A90-1 = {A90_1YR_PROPIO:.4f}")
print("  Tasas R [ev/(kg dia)] de las 4 rectas ideales:")
for (tag, k), R in R_ideal.items():
    print(f"    {tag} N_e>=k={k}: R={R:.6f}")
print(f"  A90 Ar+fondo (conocida/libre) a 62 kg dia: {1+np.sqrt(DCHI2_90/(S2_rate_conocida*62)):.4f} / "
      f"{1+np.sqrt(DCHI2_90/(S2_rate_libre*62)):.4f}")
print(f"  A90 Ar+fondo (conocida/libre) a 192 kg dia: {1+np.sqrt(DCHI2_90/(S2_rate_conocida*192)):.4f} / "
      f"{1+np.sqrt(DCHI2_90/(S2_rate_libre*192)):.4f}")
print(f"  -> {out}\n  -> {dest}  (md5 {md5(out)})")
