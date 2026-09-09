#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matriz de angulos ciegos NSI entre TODOS los pares de nucleos CEvNS reales
(idea 1 del usuario), mas la grafica estructural f(r) = (Z+2N)/(2Z+N) vs.
r = N/Z con cada nucleo marcado (idea 2, fusionada aqui por ser casi gratis
una vez que se tiene f(r)).

No depende de ninguna salida Fortran: solo usa (Z,N) tabulados. Implementa
en codigo lo que en doc/metodologia.tex esta derivado a mano (Sec. "El limite
es estructural", Ec. eq:blind_dir, eq:theta_blind, eq:ab_ratio, eq:fprime,
Tabla tab:estructural) -- no se re-deriva nada, solo se generaliza a una
matriz N x N y se agrega Cs (Z=55,N=78) para completar CsI de COHERENT
(la Tabla 2 de la tesis solo tenia I, no Cs).

Direccion ciega del plano NSI diagonal (eps_ee^dV, eps_ee^uV):
    d(Z,N) = (2Z+N, -(Z+2N))                              [eq:blind_dir]
    theta(Z,N) = arctan[(Z+2N)/(2Z+N)]                     [eq:theta_blind]
    theta_12 = |theta(Z1,N1) - theta(Z2,N2)|                [eq:angle]

Salidas (datos/):
    blind_angle_matrix.dat        matriz de angulos, texto
    fig_angulos_matriz.png        heatmap N x N
    fig_angulos_fr_vs_r.png       f(r) vs r=N/Z con nucleos marcados
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9, "axes.grid": True,
    "grid.color": "0.82", "grid.linewidth": 0.5,
    "xtick.direction": "in", "ytick.direction": "in",
    "legend.frameon": False, "legend.fontsize": 9,
    "lines.linewidth": 1.7,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})

# --------------------------------------------------------------- nucleos
# (nombre, Z, N). Z,N de Xe/Ar/Ge tomados directo de los constants.f90 del
# repo (fuente de verdad de la simulacion); O/Na/I de la Tabla 2 ya escrita
# en metodologia.tex; Cs NUEVO (isotopo estable unico, Z=55,N=78) para
# completar el blanco CsI[Na] de COHERENT (antes solo estaba I por separado).
NUCLEI = [
    ("O-16",   8,  8),   # FORTRAN90/N_EventosCEvNS_NSI*/ no lo usa; Tabla 2 tesis
    ("Na-23", 11, 12),   # idem; tambien componente de CsI[Na] (COHERENT)
    ("Ar-40", 18, 22),   # FORTRAN90/N_EventosCEvNS_NSIAr/constants.f90:46-47
    ("Ge-73", 32, 41),   # FORTRAN90/N_EventosCEvNS_NSI/constants.f90:25-26 (CONUS+)
    ("I-127", 53, 74),   # Tabla 2 tesis; componente de CsI (COHERENT)
    ("Cs-133", 55, 78),  # NUEVO: componente de CsI (COHERENT), no estaba en la tabla
    ("Xe-131", 54, 77),  # FORTRAN90/N_EventosCEvNS_NSIXe/constants.f90:33-34
]

def theta_deg(Z, N):
    """theta(Z,N) = arctan[(Z+2N)/(2Z+N)]  (Ec. eq:theta_blind), en grados."""
    return np.degrees(np.arctan2(Z + 2.0 * N, 2.0 * Z + N))

def f_ratio(r):
    """f(r) = (1+2r)/(2+r) ,  r = N/Z   (Ec. eq:ab_ratio)."""
    return (1.0 + 2.0 * r) / (2.0 + r)

# --------------------------------------------------------------- calculo
names = [n for n, Z, N in NUCLEI]
Zs = np.array([Z for n, Z, N in NUCLEI], dtype=float)
Ns = np.array([N for n, Z, N in NUCLEI], dtype=float)
rs = Ns / Zs
thetas = theta_deg(Zs, Ns)

order = np.argsort(rs)
names_o = [names[i] for i in order]
Zs_o, Ns_o, rs_o, th_o = Zs[order], Ns[order], rs[order], thetas[order]

n = len(names_o)
ang_matrix = np.abs(np.subtract.outer(th_o, th_o))

# --------------------------------------------------------------- tabla .dat
with open(f"{BASE}/blind_angle_matrix.dat", "w", encoding="utf-8") as f:
    f.write("# Matriz de angulos ciegos NSI |theta_i - theta_j| [grados]\n")
    f.write("# theta(Z,N) = arctan[(Z+2N)/(2Z+N)]  (Ec. eq:theta_blind, metodologia.tex)\n")
    f.write("# nucleo  Z  N  r=N/Z  theta[deg]\n")
    for nm, Z, N, r, th in zip(names_o, Zs_o, Ns_o, rs_o, th_o):
        f.write(f"# {nm:8s} {Z:5.1f} {N:5.1f} {r:7.4f} {th:8.3f}\n")
    f.write("#\n")
    header = "        " + "".join(f"{nm:>10s}" for nm in names_o) + "\n"
    f.write(header)
    for i, nm in enumerate(names_o):
        row = f"{nm:8s}" + "".join(f"{ang_matrix[i, j]:10.3f}" for j in range(n))
        f.write(row + "\n")

print("=" * 78)
print(" MATRIZ DE ANGULOS CIEGOS NSI  |theta_i - theta_j|  [grados]")
print(" theta(Z,N) = arctan[(Z+2N)/(2Z+N)]   (direccion ciega, plano eps_ee^dV-eps_ee^uV)")
print("=" * 78)
print(f"{'nucleo':8s} {'Z':>5s} {'N':>5s} {'r=N/Z':>8s} {'theta':>9s}")
for nm, Z, N, r, th in zip(names_o, Zs_o, Ns_o, rs_o, th_o):
    print(f"{nm:8s} {Z:5.0f} {N:5.0f} {r:8.4f} {th:8.3f}°")
print("-" * 78)
print("Matriz |theta_i - theta_j| [grados]:")
hdr = "        " + "".join(f"{nm:>9s}" for nm in names_o)
print(hdr)
for i, nm in enumerate(names_o):
    row = f"{nm:8s}" + "".join(f"{ang_matrix[i, j]:9.2f}" for j in range(n))
    print(row)
print("-" * 78)
# pares de interes citados en la tesis, como chequeo de consistencia
def ang(a, b):
    return ang_matrix[names_o.index(a), names_o.index(b)]
print(f"Xe-Ar = {ang('Xe-131','Ar-40'):.2f}°   "
      f"Ge-Xe = {ang('Ge-73','Xe-131'):.2f}°   "
      f"Ge-Ar = {ang('Ge-73','Ar-40'):.2f}°   "
      f"O-Xe (maximo realista) = {ang('O-16','Xe-131'):.2f}°")
i_max = np.unravel_index(np.argmax(ang_matrix), ang_matrix.shape)
print(f"Par con MAYOR angulo de toda la matriz: {names_o[i_max[0]]}-{names_o[i_max[1]]} "
      f"= {ang_matrix[i_max]:.2f}°")
print(f"Rango total de theta: [{th_o.min():.2f}°, {th_o.max():.2f}°]  "
      f"(franja de {th_o.max()-th_o.min():.2f}°)")

# --------------------------------------------------------------- fig 1: heatmap
fig, ax = plt.subplots(figsize=(6.4, 5.6))
im = ax.imshow(ang_matrix, cmap="Greys", vmin=0, vmax=ang_matrix.max())
ax.set_xticks(range(n)); ax.set_xticklabels(names_o, rotation=45, ha="right")
ax.set_yticks(range(n)); ax.set_yticklabels(names_o)
for i in range(n):
    for j in range(n):
        val = ang_matrix[i, j]
        color = "white" if val > 0.6 * ang_matrix.max() else "black"
        ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color=color)
ax.set_title(r"Matriz de ángulos ciegos NSI $|\theta_i-\theta_j|$ [grados]"
             "\n(orden creciente en $N/Z$)")
cbar = fig.colorbar(im, ax=ax, shrink=0.85)
cbar.set_label(r"$|\theta_i-\theta_j|$ [grados]")
fig.tight_layout()
fig.savefig(f"{BASE}/fig_angulos_matriz.png")
print(f"\n  {BASE}/fig_angulos_matriz.png")

# --------------------------------------------------------------- fig 2: f(r) vs r
r_curve = np.linspace(1.0, 1.65, 400)
f_curve = f_ratio(r_curve)

fig2, ax2 = plt.subplots(figsize=(6.4, 5.0))
ax2.plot(r_curve, f_curve, "-", color="0.25", lw=1.8,
         label=r"$f(r)=\dfrac{Z+2N}{2Z+N}=\dfrac{1+2r}{2+r}$")
ax2.scatter(rs_o, f_ratio(rs_o), s=55, color="black", zorder=5)
# offsets manuales para separar las etiquetas de I-127/Cs-133/Xe-131, que caen
# casi encimadas (r entre 1.40 y 1.43 -- es justamente el punto de la Fig.)
offsets = {
    "O-16": (6, 6), "Na-23": (6, 6), "Ar-40": (6, 6), "Ge-73": (6, 6),
    "I-127": (-55, 8), "Cs-133": (6, -16), "Xe-131": (6, 10),
}
for nm, r, th in zip(names_o, rs_o, th_o):
    ax2.annotate(nm, (r, f_ratio(r)), textcoords="offset points",
                 xytext=offsets.get(nm, (6, 6)), fontsize=9)
ax2.axhline(1.0, color="0.7", ls=":", lw=1)
ax2.set_xlabel(r"$r = N/Z$")
ax2.set_ylabel(r"$f(r)$  (pendiente de la dirección ciega)")
ax2.set_title("Por qué ningún par de blancos CE$\\nu$NS reales rompe\n"
              "la degeneración NSI: $f(r)$ es monótona y acotada")
ax2.legend(loc="upper left")
ax2.text(0.98, 0.04,
         r"$f(1)=1,\ \ f(\infty)\to 2$" "\n"
         r"$f'(r)=\dfrac{3}{(2+r)^2}>0$ (monótona)",
         transform=ax2.transAxes, ha="right", va="bottom", fontsize=8.5,
         bbox=dict(boxstyle="round", fc="white", ec="0.6"))
fig2.tight_layout()
fig2.savefig(f"{BASE}/fig_angulos_fr_vs_r.png")
print(f"  {BASE}/fig_angulos_fr_vs_r.png")
print(f"  {BASE}/blind_angle_matrix.dat")
