#!/usr/bin/env python3
"""
Tabla de LArNEST para el argon (retroceso nuclear, 218 V/cm), con el mismo
formato que la del xenon (nest.py):
    T_nr[keV]   Qy[e-/keV]   F = Var(N_e)/<N_e>
Salida: datos/nest_Ar_218V_dense.txt, que lee Tnr_to_e.f90 (carpeta Ar).

Qy sale de lar.get_nr_yields y F se estima muestreando
lar.get_yield_fluctuations. LArNEST coincide con los cinco puntos de ReD 2025
(arXiv:2510.16404) dentro de su error (-0,06 a -1,56 sigma).

F_MODE: 'nest' (el de la tesis), 'poisson', 'xe_like' o 'const:<v>'.
RED_ANCHOR = True reescala Qy para que pase por los puntos de ReD.
"""
from pathlib import Path
import numpy as np
import nestpy

# Semilla fija para que la tabla sea reproducible
nestpy.RandomGen.rndm().set_seed(20260911)
nestpy.RandomGen.rndm().lock_seed()

# --------------------------- configuracion -------------------------------
DRIFT_V = 218.0       # V/cm  (RED-100; ReD midio a 200, diferencia < 1% en Qy)
DENSITY = 1.40        # g/cm3, LAr ~87 K
T_MIN, T_MAX, N_T = 0.15, 6.0, 4000   # malla densa para Qy (NR)
F_MODE  = 'nest'      # 'nest' | 'poisson' | 'xe_like' | 'const:<valor>'
N_FGRID = 220         # nodos donde se muestrea F (luego se interpola)
N_SAMP  = 8000        # muestras de get_yield_fluctuations por nodo
RED_ANCHOR = False    # LArNEST sin reescalar, como NEST en el xenon

# ReD 2025 (arXiv:2510.16404) Tabla 1: Qy de NR de Ar, medido, 200 V/cm.
RED_E  = np.array([2.40, 3.53, 4.52, 5.48, 7.63])          # keV (E_r medio)
RED_QY = np.array([7.42, 6.99, 6.24, 5.77, 5.08])          # e-/keV
RED_DQY = np.array([0.42, 0.34, 0.30, 0.26, 0.22])         # incert. (stat (+) sist en cuadratura)

DATOS = Path(__file__).resolve().parents[1] / "datos"
XE_TABLE  = str(DATOS / "nest_218V_dense.txt")
OUTPUT    = str(DATOS / "nest_Ar_218V_dense.txt")
# -----------------------------------------------------------------------

print(f"nestpy / NEST: {nestpy.__nest_version__}")

det = nestpy.detectors.VDetector()      # detector base: el yield no depende del detector
lar = nestpy.LArNEST(det)
NR  = nestpy.LArInteraction.NR

p = lar.get_nr_yields_parameters()
print(f"Modelo NR LArNEST: alpha={p.alpha} beta={p.beta} gamma={p.gamma} "
      f"delta={p.delta} epsilon={p.epsilon} eta={p.eta} zeta={p.zeta}")

# --------------------------- Qy(T) NR (malla densa) --------------------
T_vals   = np.linspace(T_MIN, T_MAX, N_T)
Qy_model = np.array([lar.get_nr_yields(float(T), DRIFT_V, DENSITY).Ne / T for T in T_vals])

if RED_ANCHOR:
    # correccion multiplicativa: cociente ReD/LArNEST en los 5 puntos medidos,
    # interpolado en log(E), plano fuera del rango [2.4, 7.6] keV.
    qy_at_red = np.array([lar.get_nr_yields(float(E), DRIFT_V, DENSITY).Ne / E
                          for E in RED_E])
    ratio_red = RED_QY / qy_at_red
    corr = np.interp(np.log(T_vals), np.log(RED_E), ratio_red,
                     left=ratio_red[0], right=ratio_red[-1])
    Qy_vals = Qy_model * corr
    print(f"RED_ANCHOR: cociente ReD/LArNEST en los 5 puntos = "
          f"{np.round(ratio_red, 3)}  (media {ratio_red.mean():.3f})")
    qy_note = ("Qy(NR) = LArNEST reescalado para pasar por ReD 2025 Tabla 1 "
               "(2.4-7.6 keV). T_nr<2 keV: EXTRAPOLACION DE MODELO (sin medida).")
else:
    Qy_vals = Qy_model
    qy_note = "Qy(NR) = LArNEST puro (sin anclar a ReD)"

# ------------------------------- F(T) ---------------------------------
if F_MODE == 'nest':
    Tg = np.linspace(T_MIN, T_MAX, N_FGRID)
    Fg = np.ones(N_FGRID)
    print(f"Muestreando F(T) de LArNEST: {N_FGRID} nodos x {N_SAMP} muestras...")
    for k, T in enumerate(Tg):
        y = lar.get_nr_yields(float(T), DRIFT_V, DENSITY)
        d = np.fromiter(
            (lar.get_yield_fluctuations(NR, y, DENSITY).NeFluctuation for _ in range(N_SAMP)),
            dtype=float, count=N_SAMP)
        m = d.mean()
        if m > 1e-6:
            Fg[k] = d.var() / m
    F_vals = np.interp(T_vals, Tg, Fg)
    f_note = (f"F muestreada de LArNEST get_yield_fluctuations "
              f"({N_FGRID}x{N_SAMP}); sub-Poissoniana para NR")
elif F_MODE == 'poisson':
    F_vals = np.ones(N_T)
    f_note = "F = 1.0 (Poisson); supuesto 'sin fluctuacion' de ReD"
elif F_MODE == 'xe_like':
    xt = np.loadtxt(XE_TABLE, comments='#')
    F_vals = np.interp(T_vals, xt[:, 0], xt[:, 2])
    f_note = f"F copiada de la tabla de Xe ({XE_TABLE}) - sistematico"
elif F_MODE.startswith('const:'):
    c = float(F_MODE.split(':', 1)[1])
    F_vals = np.full(N_T, c)
    f_note = f"F = {c} constante - sistematico"
else:
    raise SystemExit(f"F_MODE no reconocido: {F_MODE}")

# ------------------------------ escribir NR ---------------------------
with open(OUTPUT, "w") as f:
    f.write("# Tabla NEST para RED-100 con argon (retroceso nuclear)\n")
    f.write(f"# Modelo: nestpy.LArNEST (NEST {nestpy.__nest_version__}), "
            f"E_drift={DRIFT_V} V/cm, rho={DENSITY} g/cm3\n")
    f.write(f"# {qy_note}\n")
    f.write(f"# F_MODE = {F_MODE} : {f_note}\n")
    f.write("# T_nr[keV]   Qy[e-/keV]   F=Var(Ne)/<Ne>\n")
    f.write(f"# N_puntos = {N_T}\n")
    for T, Qy, F in zip(T_vals, Qy_vals, F_vals):
        f.write(f"{T:.6f}   {Qy:.6f}   {F:.6f}\n")

print(f"\n'{OUTPUT}' creado ({N_T} puntos). Verificacion NR:")
print("  T[keV]   Qy[e-/keV]   <Ne>     F        [<2 keV = extrapolacion]")
for Tt in [0.20, 0.25, 0.35, 0.50, 0.70, 1.0, 2.0, 3.5, 5.0]:
    j = int(np.argmin(np.abs(T_vals - Tt)))
    flag = "  (medido)" if Tt >= 2.0 else ""
    print(f"  {T_vals[j]:6.3f}   {Qy_vals[j]:9.3f}   {Qy_vals[j]*T_vals[j]:6.2f}   "
          f"{F_vals[j]:.4f}{flag}")
