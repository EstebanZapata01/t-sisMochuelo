#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cascada (waterfall) de por que el argon le gana al xenon.

La Tabla de cargas debiles dice que Xe tiene ventaja de coherencia
(x3.8 por unidad de masa); la comparacion ideal dice que Ar termina
ganando por hasta x254 con umbral N_e>=4. Este script muestra en que
paso exacto de la cadena se invierte la ventaja, panel por panel, todo
con datos que el pipeline ideal ya escribe (ningun calculo nuevo):

  (1) dR/dT           -- espectro de retroceso nuclear (teoria pura,
                          antes de NEST/LArNEST y de la extraccion).
                          Aqui manda la coherencia: Xe por encima.
  (2) R_k "creados"   -- tras el modelo binomial de Fano (yield +
                          fluctuacion de ionizacion). Xe cae en picado
                          hacia N_e alto porque su yield en el umbral
                          es ~0 y sus retrocesos son blandos.
  (3) R_k "extraidos" -- tras la extraccion de la interfase (EEE).
                          Xe (EEE ~ 0.33 efectiva sobre p_F) pierde casi
                          todo; Ar (EEE ~ 0.99) no pierde nada.
  (4) R_tot(N_e >= k) -- integrado. El punto final: la ventaja Ar/Xe.

Entradas (datos/):
  espectro_continuo{Xe,Ar}.dat   : T_nr[keV]  Kop  Mue  Comb
  ionization_electrones{,_Ar}.dat : N_e  creados  extraidos  [ev/(kg dia)]
  sensib_ideal_{Xe,Ar}.dat        : NE_LO F_mode R_tot ... (fila Fnest)

Salidas: datos/fig_waterfall_{1_retroceso,2_extraidos,3_integrado}.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
C_XE, C_AR = "#33546e", "#a86a43"



def load(path):
    return np.loadtxt(path, comments="#")


dRdT_xe = load(f"{BASE}/espectro_continuoXe.dat")   # T Kop Mue Comb
dRdT_ar = load(f"{BASE}/espectro_continuoAr.dat")
ne_xe = load(f"{BASE}/ionization_electrones.dat")    # N_e creados extraidos
ne_ar = load(f"{BASE}/ionization_electrones_Ar.dat")


def rtot_ge(path):
    """R_tot(N_e>=k) para k=1..4 de la fila Fnest de sensib_ideal_*."""
    k, r = [], []
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            if p[1] == "Fnest":
                k.append(int(p[0]))
                r.append(float(p[2]))
    return np.array(k), np.array(r)


k_xe, rtot_xe = rtot_ge(f"{BASE}/sensib_ideal_Xe.dat")
k_ar, rtot_ar = rtot_ge(f"{BASE}/sensib_ideal_Ar.dat")

from estilo_tesis import FIG15
aplicar(grande=True)

def nueva():
    return plt.subplots(figsize=FIG15)

def guarda(fig, nombre):
    fig.tight_layout(); fig.savefig(f"{BASE}/{nombre}.png"); plt.close(fig)
    print(f"  -> {BASE}/{nombre}.png")

# ---- (1): dR/dT, teoria pura, con la ROI sombreada
fig, a = nueva()
a.plot(dRdT_xe[:, 0], dRdT_xe[:, 3], color=C_XE, label="Xe")
a.plot(dRdT_ar[:, 0], dRdT_ar[:, 3], color=C_AR, label="Ar")
a.axvspan(0.2, 1.0, color="0.9", zorder=0)
a.set_yscale("log"); a.set_xlim(0, 4)
a.set_xlabel(r"$T_{\rm nr}$ [keV]")
a.set_ylabel(r"$dR/dT$ [ev/(kg$\cdot$día$\cdot$keV)]")
a.legend(loc="upper right")
a.text(0.6, 0.06, "ROI CE$\\nu$NS\n(0,2–1 keV)", transform=a.get_xaxis_transform(),
       ha="center", va="bottom", fontsize=10.5, color="0.35")
guarda(fig, "fig_waterfall_1_retroceso")

# ---- (2): extraidos (creados punteados)
fig, a = nueva()
a.step(ne_xe[:, 0], ne_xe[:, 2], where="mid", color=C_XE, label="Xe, extraídos")
a.step(ne_ar[:, 0], ne_ar[:, 2], where="mid", color=C_AR, label="Ar, extraídos")
a.step(ne_xe[:, 0], ne_xe[:, 1], where="mid", color=C_XE, lw=1.2, ls=":", alpha=0.7, label="Xe, creados")
a.step(ne_ar[:, 0], ne_ar[:, 1], where="mid", color=C_AR, lw=1.2, ls=":", alpha=0.7, label="Ar, creados")
a.set_yscale("log"); a.set_xlim(0, 15); a.set_ylim(1e-4, 5e1)
a.set_xlabel(r"$N_e$")
a.set_ylabel(r"$R_k$ [ev/(kg$\cdot$día)]")
a.legend(loc="upper right")
guarda(fig, "fig_waterfall_2_extraidos")

# ---- (3): R_tot(N_e>=k) integrado + cociente Ar/Xe
fig, a = nueva()
a.plot(k_xe, rtot_xe, "o-", color=C_XE, label="Xe")
a.plot(k_ar, rtot_ar, "o-", color=C_AR, label="Ar")
a.set_yscale("log"); a.set_xticks([1, 2, 3, 4]); a.set_xlim(0.7, 4.5)
a.set_ylim(rtot_xe.min() * 0.4, rtot_ar.max() * 8)
a.set_xlabel(r"umbral $N_e \geq k$")
a.set_ylabel(r"$R_{\rm tot}$ [ev/(kg$\cdot$día)]")
a.legend(loc="lower left")
ratio = rtot_ar / rtot_xe
for kk, rr in zip(k_xe, ratio):
    a.annotate(rf"Ar/Xe $\times${rr:.0f}", (kk, rtot_ar[kk - 1]), textcoords="offset points",
               xytext=(0, 10), ha="center", fontsize=11, color="0.25")
guarda(fig, "fig_waterfall_3_integrado")

print("=" * 74)
print(" CASCADA Xe vs Ar")
print("=" * 74)
print(f"  R_tot(N_e>=1): Xe {rtot_xe[0]:.2f}   Ar {rtot_ar[0]:.2f}   Ar/Xe x{ratio[0]:.1f}")
print(f"  R_tot(N_e>=4): Xe {rtot_xe[3]:.3f}  Ar {rtot_ar[3]:.2f}   Ar/Xe x{ratio[3]:.0f}")
print(f"  Xe pierde en la extracción (N_e=4): creados {ne_xe[4,1]:.3f} -> extraídos {ne_xe[4,2]:.4f}")
print(f"  Ar conserva (N_e=4):                creados {ne_ar[4,1]:.3f} -> extraídos {ne_ar[4,2]:.3f}")
