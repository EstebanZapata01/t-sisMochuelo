#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
La forma cuadrica general de la degeneracion NSI (Ec. blind_conic) es una
CONICA. El caso por defecto ipar=5 la hace una circunferencia; otros pares
de parametros dan una elipse propia. Esto muestra que la generalidad de la
Ec. no es solo algebra: se ve.

  ipar=5  (eps_ee^dV vs eps_emu^dV):
     q_ee  = eps_x (Z+2N) ,  q_emu = eps_y (Z+2N)
     blind:  (Q_W + eps_x (Z+2N))^2 + (eps_y (Z+2N))^2 = Q_W^2   -> circunferencia
             centro (rho, 0), radio rho ;  rho = |Q_W|/(Z+2N)

  ipar=6  (eps_ee^uV vs eps_emu^dV):
     q_ee  = eps_x (2Z+N) ,  q_emu = eps_y (Z+2N)
     blind:  (Q_W + eps_x (2Z+N))^2 + (eps_y (Z+2N))^2 = Q_W^2   -> ELIPSE
             centro (a_x, 0), semiejes a_x = |Q_W|/(2Z+N), a_y = |Q_W|/(Z+2N)

El caso ipar=5 se contrasta con su contorno NUMERICO Delta_chi2 = 4.605 de
la grilla ideal (chi2_nsi_2DXe_ideal.dat): la curva analitica cae sobre el
borde de la region numerica. (El ipar=6 se muestra analitico; el driver
ideal solo barre ipar=5, y con el dato REAL de Xe toda la caja |eps|<=1
queda permitida — A_90~107 > max A_amp~74 — asi que no hay contorno cerrado
que dibujar de esos.)

Entrada: datos/chi2_nsi_2DXe_ideal.dat
Salida : datos/fig_ellipse_vs_circle.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
Z, N = 54.0, 77.0
S2W = 0.23857
QW = -N / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * Z
ZN2, Z2N = Z + 2 * N, 2 * Z + N
C5, C6 = "#33546e", "#8f4444"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.fontsize": 8.6,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})

ph = np.linspace(0, 2 * np.pi, 500)
# ipar=5: circunferencia
rho = abs(QW) / ZN2
cx5, cy5 = rho + rho * np.cos(ph), rho * np.sin(ph)
# ipar=6: elipse
ax_, ay_ = abs(QW) / Z2N, abs(QW) / ZN2
cx6, cy6 = ax_ + ax_ * np.cos(ph), ay_ * np.sin(ph)

fig, ax = plt.subplots(figsize=(7.0, 5.4))

# contorno numerico ideal (ipar=5)
d = np.loadtxt(f"{BASE}/chi2_nsi_2DXe_ideal.dat", comments="#")
x = np.unique(d[:, 0]); y = np.unique(d[:, 1])
G = d[:, 2].reshape(len(y), len(x))
ax.contour(x, y, G, levels=[4.605], colors=["0.55"], linewidths=1.1)
ax.plot([], [], color="0.55", lw=1.1, label=r"ipar=5 numérico ($\Delta\chi^2=4{,}605$, ideal)")

ax.plot(cx5, cy5, color=C5, lw=2.0, ls="--",
        label=fr"ipar=5 analítico: circunferencia ($\rho={rho:.3f}$)")
ax.plot(cx6, cy6, color=C6, lw=2.0, ls="-",
        label=fr"ipar=6 analítico: elipse ($a_x={ax_:.3f}$, $a_y={ay_:.3f}$)")
ax.plot(0, 0, "k+", ms=12, mew=1.6, label="SM ($\\varepsilon=0$)")

ax.set_aspect("equal")
ax.set_xlim(-0.15, 0.55)
ax.set_ylim(-0.28, 0.28)
ax.set_xlabel(r"$\varepsilon_x$  (ipar=5: $\varepsilon_{ee}^{dV}$;  ipar=6: $\varepsilon_{ee}^{uV}$)")
ax.set_ylabel(r"$\varepsilon_y = \varepsilon_{e\mu}^{dV}$")
ax.set_title("La degeneración es una cónica: círculo (ipar=5) o elipse (ipar=6)")
ax.legend(loc="upper right")

fig.tight_layout()
fig.savefig(f"{BASE}/fig_ellipse_vs_circle.png")

print("=" * 66)
print(f"  Q_W = {QW:.3f}   Z+2N = {ZN2:.0f}   2Z+N = {Z2N:.0f}")
print(f"  ipar=5 circunferencia: centro ({rho:.4f}, 0), radio {rho:.4f}")
print(f"  ipar=6 elipse:        centro ({ax_:.4f}, 0), semiejes "
      f"a_x={ax_:.4f}  a_y={ay_:.4f}  (excentricidad {np.sqrt(1-(ay_/ax_)**2):.3f})")
print(f"\n  {BASE}/fig_ellipse_vs_circle.png")
