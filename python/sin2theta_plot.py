#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reinterpretacion de la sensibilidad como cota sobre sin^2(theta_W).

sin^2(theta_W) entra en la prediccion CEvNS solo via la carga debil
Q_W(s) = -N/2 + (1-4s)/2 * Z, asi que la amplitud efectiva respecto al SM
es A(s) = [Q_W(s)/Q_W(s0)]^2 (s0 = 0.23857). Reparametrizando el chi2(A)
ya validado se obtiene chi2(s) sin rederivar nada.

Panel A -- Xe, DATO REAL (motor generico + generic_input_Xe.dat):
  todo s in [0, 0.5] queda con Delta_chi2 < 2.706 -> el residuo real de
  Xe NO restringe sin^2(theta_W). Ademas A(s) in [~0.095, ~10.3] para
  s in [0,1]: un corrimiento puro del angulo de Weinberg NUNCA lleva a
  A = 0, mientras que el A_best ~ 0 medido si es alcanzable por una NSI
  bien elegida.

Panel B -- Xe y Ar IDEALES (Asimov, chi2_ideal_nest, N_e>=1):
  un run ideal de Ar acota sin^2(theta_W) a ~+-1.2% del SM; el de Xe a
  ~+-2.6%. La ventaja de Ar (umbral bajo -> mas eventos) tambien aplica
  al canal electrodebil.

Panel C -- validacion del metodo con CONUS+:
  el motor generico + s-scan reproduce el programa chi2_sin2theta
  (mismo dato de CONUS+, distinta implementacion).

Entradas (datos/):
  generic_Xe_sin2theta.dat        : s  A(s)  chi2  dchi2      (Xe real)
  chi2_sin2theta_ideal_{Xe,Ar}.dat: NE_LO  s  A(s)  dchi2     (Asimov)
  generic_conus_sin2theta.dat     : s  A(s)  chi2  dchi2      (CONUS+ nuevo)
  chi2_sin2theta.dat              : s  chi2                    (CONUS+ standalone)

Salida: datos/fig_sin2theta_XeAr.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
S0 = 0.23857
CL = [(1.00, r"$1\sigma$"), (2.706, r"$90\%$"), (3.84, r"$2\sigma$")]


def cl_lines(ax):
    ymax = ax.get_ylim()[1]
    for y, lab in CL:
        if y > ymax:
            continue
        ax.axhline(y, color="0.55", lw=0.8, ls="-" if y == 2.706 else ":")
        ax.text(1.005, y, lab, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=7.6, color="0.4")


def band(dchi2, s):
    m = dchi2 <= 2.706
    return (s[m].min(), s[m].max()) if m.any() else (np.nan, np.nan)


fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.2))

# ---------- Panel A: Xe real
d = np.loadtxt(f"{BASE}/generic_Xe_sin2theta.dat", comments="#")
s, As, dchi2 = d[:, 0], d[:, 1], d[:, 3]
a = ax[0]
a.set_xlim(0, 0.5)
a.set_ylim(0, 3.2)
a.plot(s, dchi2, color=C_XE, lw=2.0)
a.axvline(S0, color="0.4", ls="--", lw=1.0)
cl_lines(a)
a.set_xlabel(r"$\sin^2\theta_W$")
a.set_ylabel(r"$\Delta\chi^2(s)$")
a.set_title("(A) Xe, dato real: sin restricción")
a.text(0.04, 0.90, r"todo $s\in[0,0.5]$ permitido",
       transform=a.transAxes, va="top", fontsize=8.4)
a.text(S0 + 0.015, 2.55, "SM", ha="left", fontsize=8, color="0.4")

# ---------- Panel B: Xe & Ar ideal (N_e>=1). Eje x = desviacion % respecto al SM
a = ax[1]
a.set_xlim(-4, 4)
a.set_ylim(0, 5)
for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
    di = np.loadtxt(f"{BASE}/chi2_sin2theta_ideal_{tag}.dat", comments="#")
    m1 = di[:, 0] == 1
    si, dci = di[m1, 1], di[m1, 3]
    lo, hi = band(dci, si)
    half = 100.0 * (hi - lo) / 2 / S0
    a.plot(100.0 * (si - S0) / S0, dci, color=col,
           label=fr"{tag} ideal:  $\pm{half:.1f}\%$")
    a.axvspan(100 * (lo - S0) / S0, 100 * (hi - S0) / S0, color=col, alpha=0.14)
a.axvline(0.0, color="0.4", ls="--", lw=1.0)
cl_lines(a)
a.set_xlabel(r"desviación de $\sin^2\theta_W$ respecto al SM  [%]")
a.set_title(r"(B) Proyección Asimov ($N_e\geq1$)")
a.legend(loc="upper center")
a.text(0.12, 0.35, "SM", ha="left", fontsize=8, color="0.4")

# ---------- Panel C: validacion CONUS+
a = ax[2]
a.set_xlim(0.13, 0.35)
a.set_ylim(0, 6)
dg = np.loadtxt(f"{BASE}/generic_conus_sin2theta.dat", comments="#")
sg, dcg = dg[:, 0], dg[:, 3]
ds_ = np.loadtxt(f"{BASE}/chi2_sin2theta.dat", comments="#")
ss, cs = ds_[:, 0], ds_[:, 1]
dcs = cs - cs.min()
a.plot(ss, dcs, color="black", lw=2.2, label="chi2_sin2theta (directo)")
a.plot(sg, dcg, color="#8f4444", lw=1.2, ls="--", label="motor genérico (s-scan)")
a.axvline(S0, color="0.4", ls="--", lw=1.0)
cl_lines(a)
a.set_xlabel(r"$\sin^2\theta_W$")
a.set_title("(C) Validación del método con CONUS+")
a.legend(loc="upper center")
sbest_g = sg[np.argmin(dcg)]
sbest_s = ss[np.argmin(dcs)]
a.text(0.97, 0.06, f"$s_{{\\rm best}} = {sbest_s:.3f}$ vs $ {sbest_g:.3f}$;  "
       f"$\\chi^2_{{\\min}}$ idéntico",
       transform=a.transAxes, va="bottom", ha="right", fontsize=7.8, color="0.35")

fig.suptitle(r"Cota sobre $\sin^2\theta_W$ vía $A(s)=[Q_W(s)/Q_W(s_0)]^2$",
             y=1.02, fontsize=12)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_sin2theta_XeAr.png")

print("=" * 74)
print(" REINTERPRETACION sin^2(theta_W)")
print("=" * 74)
print(f"  Xe real   : todo s in [{band(dchi2,s)[0]:.3f}, {band(dchi2,s)[1]:.3f}] "
      f"(= rango escaneado -> sin restriccion);  A(s) min = {As.min():.4f}")
for tag in ("Xe", "Ar"):
    di = np.loadtxt(f"{BASE}/chi2_sin2theta_ideal_{tag}.dat", comments="#")
    for nl in (1, 4):
        m = di[:, 0] == nl
        lo, hi = band(di[m, 3], di[m, 1])
        print(f"  {tag} ideal N_e>={nl}: [{lo:.4f}, {hi:.4f}]  (+-{100*(hi-lo)/2/S0:.1f}% de s0)")
print(f"  CONUS+ check: s_best directo {sbest_s:.4f}  vs  generico {sbest_g:.4f}")
print(f"\n  {BASE}/fig_sin2theta_XeAr.png")
