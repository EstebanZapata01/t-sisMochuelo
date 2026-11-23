#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nsi_core.py -- funciones de la NSI vectorial (las mismas formulas que el
select case(ipar) de los programas Fortran), para barrer tambien planos con
dos sabores de quark del mismo par de leptones.

Q_W = -N/2 + (1-4 s2w)/2 Z. Un eps_ea^{fV} suma eps (2Z+N) [uV] o
eps (Z+2N) [dV] a q_a, y A_amp = q_eff^2/Q_W^2 con
q_eff^2 = (Q_W+q_ee)^2 + q_emu^2 + q_etau^2.
"""
import numpy as np

S2W = 0.23857

TARGET = {
    "Xe": dict(Z=54.0, N=77.3879, roi=(4, 7), expo_base=192.0),
    "Ar": dict(Z=18.0, N=22.0,    roi=(1, 5), expo_base=62.0),
    # Ge natural de CONUS+ (abundancias de De Romeri et al.)
    "Ge": dict(Z=32.0, N=0.2057 * 38 + 0.2745 * 40 + 0.0775 * 41 + 0.3650 * 42 + 0.0773 * 44),
}
for t in TARGET.values():
    t["QW"] = -t["N"] / 2.0 + (1.0 - 4.0 * S2W) / 2.0 * t["Z"]

COEF = {"dV": lambda Z, N: Z + 2.0 * N, "uV": lambda Z, N: 2.0 * Z + N}


def q_lepton(eps, flavor, Z, N):
    """Contribucion de un eps_{e\\alpha}^{fV} (flavor='dV'|'uV') a q_alpha."""
    return eps * COEF[flavor](Z, N)


def A_amp(QW, q_ee=0.0, q_emu=0.0, q_etau=0.0):
    """Amplitud de tasa CEvNS relativa al SM; A_amp=1 <=> q_eff^2=Q_W^2."""
    return ((QW + q_ee) ** 2 + q_emu ** 2 + q_etau ** 2) / QW ** 2


def A_amp_grid(ex, ey, Z, N, QW, flavor_x, flavor_y, lepton_x="ee", lepton_y="emu"):
    """A_amp en una grilla (ex, ey) para dos ejes eps^{flavor}_{leptonX,Y}.
    Si lepton_x == lepton_y (p.ej. los dos son 'ee'), ambos suman al MISMO
    q_ee (plano de dos sabores de quark sobre el mismo leptón)."""
    q = {"ee": 0.0, "emu": 0.0, "etau": 0.0}
    q[lepton_x] = q[lepton_x] + q_lepton(ex, flavor_x, Z, N)
    q[lepton_y] = q[lepton_y] + q_lepton(ey, flavor_y, Z, N)
    return A_amp(QW, q["ee"], q["emu"], q["etau"])


def leer_espectro_ideal(base, tag, modo="Fnest"):
    """-> N_e (int array), R_bin [ev/(kg dia)] de datos/espectro_Ne_ideal_{tag}.dat."""
    d = np.loadtxt(f"{base}/espectro_Ne_ideal_{tag}.dat", comments="#")
    col = 1 if modo == "Fnest" else 2
    return d[:, 0].astype(int), d[:, col]


def S2_roi(base, tag, roi=None, modo="Fnest"):
    """S2 = Sum_{k en ROI} R_k  [ev/(kg dia)]; S2(t)=t*S2_roi (Asimov, sin fondo).
    roi=(k_lo,k_hi) inclusive; None -> el ROI real ya establecido en TARGET."""
    Ne, R = leer_espectro_ideal(base, tag, modo)
    lo, hi = roi or TARGET[tag]["roi"]
    return float(R[(Ne >= lo) & (Ne <= hi)].sum())


def S2_ideal(base, tag, k, modo="Fnest"):
    """S2 = Sum_{Ne>=k} R_k, sin corte superior ni fondo (como sensib_ideal_*.dat)."""
    Ne, R = leer_espectro_ideal(base, tag, modo)
    return float(R[Ne >= k].sum())


def Rk_independiente(base, tag, A_amp=1.0, kmax=15):
    """R_k [ev/(kg dia)] calculado en Python desde espectro_continuo{tag}.dat
    con el modelo binomial, para comparar bin a bin con chi2_ideal_nest.f90."""
    from scipy.stats import binom
    d = np.loadtxt(f"{base}/espectro_continuo{tag}.dat", comments="#")
    T, dRdT = d[:, 0], d[:, 3]
    dT = T[1] - T[0]
    nest = np.loadtxt(f"{base}/nest_{'Ar_' if tag == 'Ar' else ''}218V_dense.txt", comments="#")
    Qy = np.interp(T, nest[:, 0], nest[:, 1])
    F = np.interp(T, nest[:, 0], nest[:, 2])
    EEE = 0.328 if tag == "Xe" else 0.99
    Rk = np.zeros(kmax + 1)
    for i in range(len(T)):
        if dRdT[i] <= 0:
            continue
        lam = Qy[i] * T[i]
        pF = 1.0 - F[i]
        if lam <= 0 or pF <= 1e-6:
            continue
        nF = int(round(lam / pF))
        ks = np.arange(0, min(nF, kmax) + 1)
        pmf = binom.pmf(ks, nF, pF * EEE)
        Rk[ks] += A_amp * dRdT[i] * dT * pmf
    return Rk   # indice = N_e verdadero (0..kmax)


def punto_en_circulo(A_target, QW, Z, N, angulo_deg, flavor="dV", lepton_x="ee", lepton_y="emu"):
    """Punto (ex, ey) de la circunferencia A_amp = A_target en la direccion
    angulo_deg del plano (q_x, q_y)."""
    a = np.radians(angulo_deg)
    qx = QW * (np.sqrt(A_target) * np.cos(a) - 1.0)
    qy = QW * np.sqrt(A_target) * np.sin(a)
    cx = COEF[flavor](Z, N)
    return qx / cx, qy / cx


# --------------------------------------------------- respuesta instrumental de Xe
SEG, SIG1 = 27.0, 7.6
PE_LO, PE_HI = 110.0, 189.0
EFF_ROI_XE = {4: 0.1369, 5: 0.3271, 6: 0.6105, 7: 0.7373}


def prob_migracion(j, k):
    from scipy.special import erf
    a = max((j - 0.5) * SEG, PE_LO); b = min((j + 0.5) * SEG, PE_HI)
    mu, sg = k * SEG, np.sqrt(k) * SIG1
    return 0.5 * (erf((b - mu) / (sg * np.sqrt(2))) - erf((a - mu) / (sg * np.sqrt(2))))


def respuesta_Xe_roi(R_true, Ne_true, roi=(4, 7)):
    """Pasa R(N_e verdadero) a R_j reconstruido en la ROI del xenon (migracion
    en PE y eficiencia), como red100PE.f90."""
    lo, hi = roi
    Rj = np.zeros(hi - lo + 1)
    for idx, j in enumerate(range(lo, hi + 1)):
        Rj[idx] = EFF_ROI_XE[j] * sum(prob_migracion(j, k) * R_true[Ne_true == k].sum()
                                        for k in Ne_true)
    return np.arange(lo, hi + 1), Rj
