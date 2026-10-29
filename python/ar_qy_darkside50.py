#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q_y(T_nr) de retroceso nuclear del argon segun el modelo de DarkSide-50
(PRD 104, 082005 (2021), arXiv:2107.08087, Sec. V, Ecs. 9-13), evaluado al
campo de deriva de la cadena (218 V/cm) en vez de los 200 V/cm de su ajuste
(el propio articulo dice que ese cambio mueve Qy < 1%, igual que LArNEST):

  Qy = (1-r) N_i / E,   1-r = ln(1+gamma N_i)/(gamma N_i),  gamma = C_box/F,
  N_i = beta * kappa(eps),  kappa(eps) = eps * s_e(eps) / (s_n(eps)+s_e(eps)),
  eps = 0.0135 * E[keV],  s_e(eps) = 0.145 * sqrt(eps)  (F(v/v0)=1, su eleccion
  central), s_n(eps) = ln(1+1.1383 fZ eps) / (2[fZ eps + 0.01321 (fZ eps)^0.21226
  + 0.19593 (fZ eps)^0.5]),  fZ = 0.953 (frenado nuclear de Ziegler).

Mejor ajuste global a 200 V/cm: C_box = 8.1 (+0.1/-0.2) V/cm,
beta = 6.8e3 (+0.1/-0.3)e3. Calibrado por sus datos AmC in-situ hasta
435 (+47/-34) eV_nr; SCENE/ARIS externos cubren ~7-118 keV_nr.

Comparacion de sustento: contra LArNEST-ReD (nest_Ar_218V_dense.txt, anclado
a ReD arXiv:2510.16404 Tabla 1, 2.4-7.6 keV) y contra los 5 puntos de ReD
directamente. Si el modelo DS-50 no es compatible (dentro de ~2-3 sigma) con
esos puntos medidos, NO se usa como sistematico (lo dice el pie de cada
figura/cifra que lo use).

Si esta sustentado, corre chi2_ideal y chi2_bkg (Ar) tres veces (central,
+1sigma, -1sigma de C_box/beta tomados como banda simple, no correlacionada)
sustituyendo solo la columna Qy de nest_Ar_218V_dense.txt (conserva la F de
LArNEST, restaura el original al salir) y compara tasas/A90 con el baseline.

Entrada: datos/nest_Ar_218V_dense.txt
Salida : datos/qy_darkside50_Ar.dat, datos/fig_qy_darkside50.png,
         stdout con las cifras para cifras_tesis_helm.md
"""
import os
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_AR, GRIS, NEGRO, FIG15, coma
aplicar(grande=True)

ROOT = os.path.abspath("..")
BASE = f"{ROOT}/datos"
DIR_AR = f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"
SCRATCH = "/tmp/ar_qy_ds50_backup"

# --------------------------------------------------------- modelo DS-50
F_FIELD = 218.0          # V/cm, campo de la cadena (su ajuste es a 200)
FZ = 0.953


def s_n_ziegler(eps):
    fe = FZ * eps
    return np.log(1.0 + 1.1383 * fe) / (2.0 * (fe + 0.01321 * fe**0.21226 + 0.19593 * fe**0.5))


def s_e_bezrukov(eps):
    return 0.145 * np.sqrt(eps)


def Qy_DS50(E_keV, C_box, beta):
    eps = 0.0135 * E_keV
    se, sn = s_e_bezrukov(eps), s_n_ziegler(eps)
    kappa = eps * se / (sn + se)
    Ni = beta * kappa
    gamma = C_box / F_FIELD
    gNi = gamma * Ni
    uno_menos_r = np.where(gNi > 1e-9, np.log1p(gNi) / np.maximum(gNi, 1e-300), 1.0)
    return uno_menos_r * Ni / E_keV


PARES = {"central": (8.1, 6.8e3), "+1sigma": (8.2, 6.9e3), "-1sigma": (7.9, 6.5e3)}

# --------------------------------------------------------- comparacion
T, Qy_base, F_base = np.loadtxt(f"{BASE}/nest_Ar_218V_dense.txt", comments="#", unpack=True)
Qy_ds50 = {k: Qy_DS50(T, *v) for k, v in PARES.items()}

RED_E = np.array([2.40, 3.53, 4.52, 5.48, 7.63])
RED_QY = np.array([7.42, 6.99, 6.24, 5.77, 5.08])
RED_DQY = np.array([0.42, 0.34, 0.30, 0.26, 0.22])
ds50_at_red = Qy_DS50(RED_E, *PARES["central"])
pull_red = (ds50_at_red - RED_QY) / RED_DQY
larnest_at_red = np.interp(RED_E, T, Qy_base)

print("=" * 78)
print(" Qy DarkSide-50 (218 V/cm) vs ReD 2025 medido y vs LArNEST-ReD (sustento)")
print("=" * 78)
for i, E in enumerate(RED_E):
    print(f"  E={E:5.2f} keV  ReD={RED_QY[i]:.2f}+-{RED_DQY[i]:.2f}  "
          f"LArNEST-ReD={larnest_at_red[i]:.3f}  DS50={ds50_at_red[i]:.3f}  "
          f"pull(DS50-ReD)={pull_red[i]:+.2f} sigma")
max_pull = np.max(np.abs(pull_red))
print(f"  max |pull| frente a ReD = {max_pull:.2f} sigma")
SUSTENTADO = max_pull < 3.0
print(f"  SUSTENTADO (<3 sigma en todos los puntos medidos): {SUSTENTADO}")

roi = (T >= 0.1) & (T <= 1.0)
dev_roi = Qy_ds50["central"][roi] / Qy_base[roi] - 1.0
print(f"\n  En la ROI CEvNS (T=0.1-1 keV): DS50/LArNEST-ReD - 1 va de "
      f"{dev_roi.min():+.1%} a {dev_roi.max():+.1%}")

# ------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=FIG15)
ax.fill_between(T, Qy_ds50["-1sigma"], Qy_ds50["+1sigma"], color=C_AR, alpha=0.18, lw=0,
                 label=r"DS-50, banda $1\sigma$ ($C_{\rm box},\beta$)")
ax.plot(T, Qy_ds50["central"], color=C_AR, lw=2.0, label="DS-50 (central)")
ax.plot(T, Qy_base, color=NEGRO, lw=1.6, ls="--", label="LArNEST-ReD (base de la cadena)")
ax.errorbar(RED_E, RED_QY, yerr=RED_DQY, fmt="o", color=NEGRO, ms=6, capsize=3,
            label="ReD 2025 (medido)")
ax.axvspan(0.1, 1.0, color=GRIS, alpha=0.15, lw=0)
ax.text(0.55, 0.04, "ROI CE$\\nu$NS", transform=ax.get_xaxis_transform(), ha="center",
        va="bottom", fontsize=10.5, color="0.35")
ax.set_xlim(0, 7.8)
ax.set_xlabel(r"$T_{\rm nr}$  [keV]")
ax.set_ylabel(r"$Q_y(T_{\rm nr})$  [e$^-$/keV]")
ax.legend(loc="lower right", fontsize=9.6)
fig.tight_layout()
out = f"{BASE}/fig_qy_darkside50.png"
fig.savefig(out)
print(f"\n  -> {out}")

with open(f"{BASE}/qy_darkside50_Ar.dat", "w", encoding="utf-8") as f:
    f.write("# T_nr[keV]  Qy_DS50_central  Qy_DS50_+1sigma  Qy_DS50_-1sigma  Qy_LArNEST_ReD\n")
    for i in range(len(T)):
        f.write(f"{T[i]:.6f}  {Qy_ds50['central'][i]:.6f}  {Qy_ds50['+1sigma'][i]:.6f}  "
                f"{Qy_ds50['-1sigma'][i]:.6f}  {Qy_base[i]:.6f}\n")

if not SUSTENTADO:
    print("\n  NO SUSTENTADO: no se corre la cadena de Ar con este Qy. Avisar a la tesis.")
    raise SystemExit(0)

# ------------------------------------------------- correr la cadena de Ar
os.makedirs(SCRATCH, exist_ok=True)
tabla_path = f"{BASE}/nest_Ar_218V_dense.txt"
shutil.copy(tabla_path, f"{SCRATCH}/nest_Ar_218V_dense_ORIGINAL.txt")

resultados = {}
try:
    for variante in ("central", "+1sigma", "-1sigma"):
        with open(tabla_path, "w", encoding="utf-8") as f:
            f.write(f"# Ar NR Qy sustituido por DarkSide-50 ({variante}), F de LArNEST\n")
            for i in range(len(T)):
                f.write(f"{T[i]:.6f}  {Qy_ds50[variante][i]:.6f}  {F_base[i]:.6f}\n")
        out_ideal = subprocess.run(["./chi2_ideal"], cwd=DIR_AR, capture_output=True, text=True)
        out_bkg = subprocess.run(["./chi2_bkg"], cwd=DIR_AR, capture_output=True, text=True)
        espectro = np.loadtxt(f"{BASE}/espectro_Ne_ideal_Ar.dat", comments="#")
        Rk = {int(row[0]): row[1] for row in espectro}   # col 1 = R_bin_Fnest (F de NEST)
        resultados[variante] = dict(Rk=Rk, ideal=out_ideal.stdout, bkg=out_bkg.stdout)
        print(f"\n  -- variante {variante} (C_box={PARES[variante][0]}, beta={PARES[variante][1]:.1f}) --")
        rate_k = {k: sum(v for ne, v in Rk.items() if ne >= k) for k in (1, 2, 3, 4)}
        print(f"    R_tot(N_e>=k) [ev/(kg dia)]: " +
              ", ".join(f"k={k}:{rate_k[k]:.3f}" for k in (1, 2, 3, 4)))
        for ln in out_ideal.stdout.splitlines():
            s = ln.strip()
            if s.startswith("---") or s.startswith("NE_LO") or (s[:1].isdigit() and len(s.split()) == 5):
                print("    [ideal] " + s)
        for ln in out_bkg.stdout.splitlines():
            if "A_90" in ln or "A90" in ln:
                print("    [bkg]   " + ln.strip())
finally:
    shutil.copy(f"{SCRATCH}/nest_Ar_218V_dense_ORIGINAL.txt", tabla_path)
    # IMPORTANTE: restaurar el .txt no alcanza -- sensib_ideal_Ar.dat,
    # bkg_bins_Ar.dat, etc. siguen con la salida de la ULTIMA variante
    # (-1sigma) hasta que chi2_ideal/chi2_bkg se vuelven a correr.
    subprocess.run(["./chi2_ideal"], cwd=DIR_AR, capture_output=True)
    subprocess.run(["./chi2_bkg"], cwd=DIR_AR, capture_output=True)
    print(f"\n  nest_Ar_218V_dense.txt restaurado al original (LArNEST-ReD);"
          f" chi2_ideal y chi2_bkg re-corridos para regenerar sensib_ideal_Ar.dat"
          f" y bkg_bins_Ar.dat limpios.")

with open(f"{BASE}/rtot_darkside50_Ar.dat", "w", encoding="utf-8") as f:
    f.write("# variante  R_tot(Ne>=1)  R_tot(Ne>=2)  R_tot(Ne>=3)  R_tot(Ne>=4)  [ev/(kg dia), F nest]\n")
    for variante in ("central", "+1sigma", "-1sigma"):
        Rk = resultados[variante]["Rk"]
        rate_k = [sum(v for ne, v in Rk.items() if ne >= k) for k in (1, 2, 3, 4)]
        f.write(f"{variante:8s} " + " ".join(f"{r:.6f}" for r in rate_k) + "\n")

print("\n  -> datos/qy_darkside50_Ar.dat")
print("  -> datos/rtot_darkside50_Ar.dat")
