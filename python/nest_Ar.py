#!/usr/bin/env python3
"""
Genera la tabla NEST para RED-100 con ARGON liquido (Ar NR, 218 V/cm).

Salida: 3 columnas   ->   datos/nest_Ar_218V_dense.txt
    T_nr[keV]   Qy[e-/keV]   F = Var(N_e)/<N_e>

Mismo formato y misma metodologia que la de Xe (nest.py -> nest_218V_dense.txt);
lo lee  FORTRAN90/N_EventosCEvNS_NSIAr/Tnr_to_e.f90.
Tesis: metodologia.tex Sec. 3.2 y Sec. 9.

--------------------------------------------------------------------------
MODELO: nestpy 2.4.5 trae la clase  `nestpy.LArNEST`  con el modelo NR de
argon liquido de NEST ("Global Analysis of Argon Yields"). NO hace falta un
"detector de ejemplo" de Ar: el yield es fisica independiente del detector,
basta el detector base VDetector().
Parametros NR (get_nr_yields_parameters): alpha=11.1 beta=0.087 gamma=0.1
delta=-0.0932 epsilon=2.998 eta=2.94 zeta=0.3.

  Qy(T):  lar.get_nr_yields(T, E_drift, rho).Ne / T     (media analitica)
  F(T) :  se MUESTREA  lar.get_yield_fluctuations(NR, y, rho).NeFluctuation
          N_SAMP veces y se toma Var/media  --  EXACTAMENTE igual que nest.py
          hace con nc.GetQuanta(...).electrons para el xenon.

VALIDACION del Qy (218 V/cm) frente a ReD arXiv:2510.16404 Tabla 1:
    E[keV]  ReD(e-/keV)  LArNEST   pull
     2.40   7.42+/-0.42   7.09    -0.78 sigma
     3.53   6.99+/-0.34   6.46    -1.56 sigma
     4.52   6.24+/-0.30   6.02    -0.73 sigma
     5.48   5.77+/-0.26   5.67    -0.38 sigma
     7.63   5.08+/-0.22   5.07    -0.06 sigma
=> compatible en 2-8 keV; LArNEST extrapola solo hasta el umbral (Lindhard
   incorporado). 200 -> 218 V/cm mueve Qy < 1%.

VALIDACION de F: la resolucion sigma(Ne)/Ne del modelo LArNEST (~12% a
Ne~10, ~8% a Ne~15) reproduce la resolucion medida por ReD (12% a Ne=10,
~7% a Ne>40). El modelo LArNEST es sub-Poissoniano para NR: F ~ 0.10-0.15
en el ROI de RED-100 (T ~ 0.3-0.8 keV), sube cerca del umbral.

F_MODE (para sistematicos):
  'nest'      -> F muestreada de LArNEST         (BASELINE de la tesis)
  'poisson'   -> F = 1                           (supuesto "sin fluctuacion")
  'xe_like'   -> F copiada de la tabla de Xe
  'const:<v>' -> F = v constante

RED_ANCHOR: si True, el Qy de retroceso NUCLEAR se reescala con una correccion
multiplicativa suave para que pase por los 5 puntos MEDIDOS de ReD 2025
(arXiv:2510.16404 Tabla 1) en 2.4-7.6 keV. El cociente ReD/LArNEST se interpola
y se mantiene plano fuera de ese rango. Por debajo de ~2 keV (el ROI CEvNS de
RED-100, N_e<=5 <-> T_nr ~ 0.1-1 keV) NO HAY MEDIDA: es LArNEST reescalado =
EXTRAPOLACION DE MODELO. Se marca como tal.

Tabla ER (nueva, OUTPUT_ER): yield de RETROCESO ELECTRONICO de LArNEST, para el
fondo beta de 39Ar (el 39Ar decae por beta -> retroceso electronico). Alimenta
N_e directamente (beta -> N_e -> extraccion), SIN paso PE. El extremo de baja
energia (E_er < 0.2 keV) es tambien extrapolacion de modelo.

39Ar: este script da el yield de senal (NR) y el yield ER para el fondo. El
espectro beta y la normalizacion (actividad) viven en chi2_bkg_nest.f90.
--------------------------------------------------------------------------
"""
import numpy as np
import nestpy

# Semilla fija (misma razon que nest.py); distinta a la de Xe.
nestpy.RandomGen.rndm().set_seed(20260911)
nestpy.RandomGen.rndm().lock_seed()

# --------------------------- configuracion -------------------------------
DRIFT_V = 218.0       # V/cm  (RED-100; ReD midio a 200, diferencia < 1% en Qy)
DENSITY = 1.40        # g/cm3, LAr ~87 K
T_MIN, T_MAX, N_T = 0.15, 6.0, 4000   # malla densa para Qy (NR)
F_MODE  = 'nest'      # 'nest' | 'poisson' | 'xe_like' | 'const:<valor>'
N_FGRID = 220         # nodos donde se MUESTREA F (luego se interpola a la malla densa)
N_SAMP  = 8000        # muestras de get_yield_fluctuations por nodo
RED_ANCHOR = True     # reescalar Qy(NR) para pasar por los 5 puntos de ReD 2025

# ReD 2025 (arXiv:2510.16404) Tabla 1: Qy de NR de Ar, medido, 200 V/cm.
RED_E  = np.array([2.40, 3.53, 4.52, 5.48, 7.63])          # keV (E_r medio)
RED_QY = np.array([7.42, 6.99, 6.24, 5.77, 5.08])          # e-/keV
RED_DQY = np.array([0.42, 0.34, 0.30, 0.26, 0.22])         # incert. (stat (+) sist en cuadratura)

XE_TABLE  = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/nest_218V_dense.txt"
OUTPUT    = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/nest_Ar_218V_dense.txt"
OUTPUT_ER = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/nest_Ar_ER_218V.txt"
EER_MIN, EER_MAX, N_EER = 0.02, 30.0, 4000   # malla ER para el fondo de 39Ar
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
    f.write("# Tabla NEST para RED-100 con ARGON  (Ar, retroceso NUCLEAR - senal CEvNS)\n")
    f.write(f"# Modelo: nestpy.LArNEST (NEST {nestpy.__nest_version__}), "
            f"E_drift={DRIFT_V} V/cm, rho={DENSITY} g/cm3\n")
    f.write(f"# {qy_note}\n")
    f.write("# ReD 2025 (arXiv:2510.16404) Tabla 1 = fuente primaria del yield NR "
            "(2.4-7.6 keV, medido).\n")
    f.write("# Bondar et al. 2017 (JINST 12, C05010): solo 80/233 keV; citado por "
            "consistencia con ref.[46], NO usado.\n")
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

# --------------------- tabla ER (fondo beta de 39Ar) -----------------
Eer_vals = np.geomspace(EER_MIN, EER_MAX, N_EER)
Qy_er = np.array([lar.get_er_yields(float(E), DRIFT_V, DENSITY).Ne / E for E in Eer_vals])
# F_er: get_fano_er si esta disponible; si no, sub-Poissoniano suave ~0.3
try:
    F_er = np.array([max(lar.get_fano_er(float(E), DENSITY), 1e-3) for E in Eer_vals])
    fer_note = "F_er = get_fano_er(E, rho) de LArNEST"
except Exception:
    F_er = np.full(N_EER, 0.3)
    fer_note = "F_er = 0.3 constante (get_fano_er no disponible) - sistematico"

with open(OUTPUT_ER, "w") as f:
    f.write("# Tabla NEST para RED-100 con ARGON  (Ar, retroceso ELECTRONICO)\n")
    f.write("# Uso: fondo beta de 39Ar (39Ar -> 39K + e- + nubar).  El espectro\n")
    f.write("#      beta y la actividad viven en chi2_bkg_nest.f90.\n")
    f.write(f"# Modelo: nestpy.LArNEST get_er_yields, E_drift={DRIFT_V} V/cm, rho={DENSITY}\n")
    f.write(f"# {fer_note}\n")
    f.write("# E_er < 0.2 keV: EXTRAPOLACION DE MODELO (sin datos ER de LAr a esa energia)\n")
    f.write("# E_er[keV]   Qy_er[e-/keV]   F_er\n")
    f.write(f"# N_puntos = {N_EER}\n")
    for E, Q, Fe in zip(Eer_vals, Qy_er, F_er):
        f.write(f"{E:.6e}   {Q:.6f}   {Fe:.6f}\n")

print(f"\n'{OUTPUT_ER}' creado ({N_EER} puntos). Verificacion ER (N_e ~ E*Qy_er):")
print("  E_er[keV]  Qy_er[e-/keV]  <Ne>")
for Ee in [0.03, 0.05, 0.08, 0.12, 0.20, 0.50, 1.0, 5.0]:
    j = int(np.argmin(np.abs(Eer_vals - Ee)))
    print(f"  {Eer_vals[j]:8.4f}   {Qy_er[j]:10.2f}   {Qy_er[j]*Eer_vals[j]:6.2f}")
