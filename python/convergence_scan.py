#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estudio de convergencia numerica: A_90 (comparacion ideal) vs. resolucion
de las mallas de integracion n_T (retroceso) y n_E (neutrino).

n_T y n_E son `parameter` en chi2_ideal_nest.f90, asi que cada punto exige
recompilar. Este script hace sed sobre la linea del parameter, recompila,
corre y recoge A_90(N_e>=1) y A_90(N_e>=4) para Xe y Ar. Deja el archivo
fuente y los binarios en su estado original (n_T=800, n_E=2000).

Si A_90 tiene un plateau claro a partir de la resolucion actual, es
evidencia de que los resultados no son artefacto de discretizacion.

Salida: datos/fig_convergencia.png  (+ datos/convergencia.dat)
"""
import os
import re
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, C_AR, C_GE, C_SM, CICLO
aplicar()

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
BASE = f"{ROOT}/datos"
DIRS = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe",
        "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}
SRC = "chi2_ideal_nest.f90"
MODS = {"Xe": "constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 "
              "Tnr_to_e.f90 mod_detector.f90",
        "Ar": "constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"}
FF = "gfortran -O2 -ffree-line-length-none"
PAT = re.compile(r"integer,\s*parameter\s*::\s*n_T\s*=\s*\d+\s*,\s*n_E\s*=\s*\d+")

NT_LIST = [200, 400, 800, 1600]
NE_LIST = [500, 1000, 2000, 4000]
BASE_NT, BASE_NE = 800, 2000

C_XE, C_AR = "#33546e", "#a86a43"
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
    "axes.linewidth": 0.9, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "legend.fontsize": 8.5, "lines.linewidth": 1.7,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})


def set_mesh(nt, ne):
    for d in DIRS.values():
        p = f"{d}/{SRC}"
        txt = open(p).read()
        new = PAT.sub(f"integer,  parameter :: n_T = {nt}, n_E = {ne}", txt, count=1)
        assert new != txt, f"patron n_T/n_E no encontrado en {p}"
        open(p, "w").write(new)


def build_run(tag):
    d = DIRS[tag]
    subprocess.run(f"{FF} -o chi2_ideal {MODS[tag]} {SRC}", cwd=d, shell=True,
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run(["./chi2_ideal"], cwd=d, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_a90(tag):
    """-> {(NE_LO): A_90_x1}  (fila Fnest)."""
    out = {}
    for ln in open(f"{BASE}/sensib_ideal_{tag}.dat"):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if p[1] == "Fnest":
            out[int(p[0])] = float(p[5])
    return out


# --- backup fuentes ---
BK = {d: open(f"{d}/{SRC}").read() for d in DIRS.values()}
results = {"nT": [], "nE": []}
try:
    # barrido en n_T (n_E fijo)
    for nt in NT_LIST:
        set_mesh(nt, BASE_NE)
        row = {"n": nt}
        for tag in ("Xe", "Ar"):
            build_run(tag)
            a = read_a90(tag)
            row[f"{tag}_ne1"] = a[1]; row[f"{tag}_ne4"] = a[4]
        results["nT"].append(row)
        print(f"  n_T={nt:5d}  Xe A90(ne1,ne4)=({row['Xe_ne1']:.4f},{row['Xe_ne4']:.4f})  "
              f"Ar=({row['Ar_ne1']:.4f},{row['Ar_ne4']:.4f})")
    # barrido en n_E (n_T fijo)
    for ne in NE_LIST:
        set_mesh(BASE_NT, ne)
        row = {"n": ne}
        for tag in ("Xe", "Ar"):
            build_run(tag)
            a = read_a90(tag)
            row[f"{tag}_ne1"] = a[1]; row[f"{tag}_ne4"] = a[4]
        results["nE"].append(row)
        print(f"  n_E={ne:5d}  Xe A90(ne1,ne4)=({row['Xe_ne1']:.4f},{row['Xe_ne4']:.4f})  "
              f"Ar=({row['Ar_ne1']:.4f},{row['Ar_ne4']:.4f})")
finally:
    # --- restaurar fuentes + rebuild baseline ---
    for d, txt in BK.items():
        open(f"{d}/{SRC}", "w").write(txt)
    for tag in ("Xe", "Ar"):
        build_run(tag)
    print("  fuentes y binarios restaurados a n_T=800, n_E=2000.")

# --- figura ---
fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.2))
for a, key, xlab in ((ax[0], "nT", r"$n_T$ (malla de retroceso)"),
                     (ax[1], "nE", r"$n_E$ (malla de neutrino)")):
    n = [r["n"] for r in results[key]]
    for tag, col in (("Xe", C_XE), ("Ar", C_AR)):
        a.plot(n, [r[f"{tag}_ne4"] for r in results[key]], "o-", color=col,
               label=fr"{tag}, $N_e\geq4$")
        a.plot(n, [r[f"{tag}_ne1"] for r in results[key]], "s--", color=col,
               alpha=0.55, label=fr"{tag}, $N_e\geq1$")
    a.axvline(BASE_NT if key == "nT" else BASE_NE, color="0.5", lw=1.0, ls=":")
    a.set_xscale("log", base=2)
    a.set_xlabel(xlab)
    a.set_ylabel(r"$A_{90}$ (ideal, $\times$SM)")
    a.grid(True, ls=":", lw=0.4, color="0.8")
a.legend = ax[0].legend(loc="center right", fontsize=8)
ax[0].set_title("Convergencia en $n_T$ ($n_E=2000$)")
ax[1].set_title("Convergencia en $n_E$ ($n_T=800$)")
fig.suptitle("Convergencia numérica: $A_{90}$ estable a partir de la "
             "resolución usada ($n_T=800$, $n_E=2000$)", y=1.02, fontsize=11.5)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_convergencia.png")

# tabla
with open(f"{BASE}/convergencia.dat", "w") as f:
    f.write("# barrido n_T (n_E=2000):  n  Xe_ne1 Xe_ne4 Ar_ne1 Ar_ne4\n")
    for r in results["nT"]:
        f.write(f"{r['n']:6d} {r['Xe_ne1']:.5f} {r['Xe_ne4']:.5f} "
                f"{r['Ar_ne1']:.5f} {r['Ar_ne4']:.5f}\n")
    f.write("# barrido n_E (n_T=800):  n  Xe_ne1 Xe_ne4 Ar_ne1 Ar_ne4\n")
    for r in results["nE"]:
        f.write(f"{r['n']:6d} {r['Xe_ne1']:.5f} {r['Xe_ne4']:.5f} "
                f"{r['Ar_ne1']:.5f} {r['Ar_ne4']:.5f}\n")

# reporte de estabilidad relativa entre la ultima resolucion y la base
def rel(seq, nbase):
    d = {r["n"]: r for r in seq}
    b = d[nbase]
    top = max(d)
    return {k: abs(d[top][k] - b[k]) / b[k] for k in ("Xe_ne4", "Ar_ne4")}

print("\n  variacion relativa A_90(N_e>=4) entre resolucion base y la mas fina:")
print(f"    n_T {BASE_NT}->{NT_LIST[-1]}: Xe {rel(results['nT'],BASE_NT)['Xe_ne4']*100:.2f}%  "
      f"Ar {rel(results['nT'],BASE_NT)['Ar_ne4']*100:.2f}%")
print(f"    n_E {BASE_NE}->{NE_LIST[-1]}: Xe {rel(results['nE'],BASE_NE)['Xe_ne4']*100:.2f}%  "
      f"Ar {rel(results['nE'],BASE_NE)['Ar_ne4']*100:.2f}%")
print(f"\n  {BASE}/fig_convergencia.png")
