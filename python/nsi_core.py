#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nsi_core.py -- piezas reutilizables de la NSI vectorial (sola tasa) para las
figuras de cierre y de combinacion Xe+Ar. Formulas identicas (mismos
coeficientes Z, 2Z+N, Z+2N) a las 15 ramas 'select case(ipar)' de
FORTRAN90/chi2_nsi_generic/chi2_nsi_generic.f90; se reimplementan aqui en
Python solo para poder barrer planos que ese select case no cubre
(dos sabores de quark sobre el MISMO par de leptones, p.ej. eps_ee^dV y
eps_ee^uV a la vez).

Un blanco (Z, N) fija Q_W = -N/2 + (1-4 s2w)/2 * Z (constants.f90 de cada
carpeta). Un eps_ee^{fV} (f=u,d) suma eps*(2Z+N) [uV] o eps*(Z+2N) [dV] a
q_ee; analogo para eps_emu, eps_etau. La amplitud de tasa (la misma para
todo N_e, ver xsections_nest.f90) es A_amp = q_eff^2/Q_W^2 con
q_eff^2 = (Q_W+q_ee)^2 + q_emu^2 + q_etau^2.
"""
import numpy as np

S2W = 0.23857

TARGET = {
    "Xe": dict(Z=54.0, N=77.3879, roi=(4, 7), expo_base=192.0),
    "Ar": dict(Z=18.0, N=22.0,    roi=(1, 5), expo_base=62.0),
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
    """S2 = Sum_{Ne>=k} R_k, SIN corte superior -- el mismo marco que
    tab:ideal (sensib_ideal_{Xe,Ar}.dat, A_90(x1)): umbral N_e>=k, eps_ROI=1,
    respuesta central, sin fondo."""
    Ne, R = leer_espectro_ideal(base, tag, modo)
    return float(R[Ne >= k].sum())


def Rk_independiente(base, tag, A_amp=1.0, kmax=15):
    """R_k [ev/(kg dia)] recalculado DESDE CERO en Python: lee
    espectro_continuo{tag}.dat (dR/dT hibrido, de red100_nest.f90 -- un
    programa Fortran DISTINTO del que genera espectro_Ne_ideal) y convoluciona
    con el modelo binomial de Fano (tabla NEST propia, interpolada), escalando
    dR/dT por A_amp (T-independiente, como exige el modelo NSI vectorial).
    Es una segunda implementacion independiente (otro archivo fuente, otra
    integracion) del mismo R_k que produce chi2_ideal_nest.f90; sirve como
    cierre bin a bin, NO se usa en ningun resultado de la tesis."""
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
    """Un punto (ex,ey) analitico sobre A_amp(ex,ey)=A_target (q_etau=0), en la
    direccion 'angulo_deg' del plano (q_x,q_y) medida desde el eje +q_x:
    q_x = QW*(sqrt(A_target)*cos(a) - 1), q_y = QW*sqrt(A_target)*sin(a)
    (verifica (QW+q_x)^2+q_y^2 = A_target*QW^2 por construccion), convertido
    a (ex,ey) dividiendo por el coeficiente de cada sabor."""
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
    """R_true(Ne_true) [ev/(kg dia)], espacio N_e VERDADERO -> R_j reconstruido
    en la ROI real de Xe (migracion en PE + eff_ROI), misma cadena lineal que
    red100PE.f90 / mod_detector.f90 (independiente de si R_true trae NSI o no:
    es una transformacion LINEAL en k, por eso preserva la factorizacion)."""
    lo, hi = roi
    Rj = np.zeros(hi - lo + 1)
    for idx, j in enumerate(range(lo, hi + 1)):
        Rj[idx] = EFF_ROI_XE[j] * sum(prob_migracion(j, k) * R_true[Ne_true == k].sum()
                                        for k in Ne_true)
    return np.arange(lo, hi + 1), Rj
