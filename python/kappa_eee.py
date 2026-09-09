#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kappa_eee.py -- factor de amplificacion logaritmica de la incertidumbre de la
eficiencia de extraccion de electrones (EEE) hasta la senal aceptada en el ROI.

    kappa_EEE = (EEE / S) * dS/dEEE ,
    S(EEE) = INT dT  w(T) * Pr[ Bin(n_F(T), p_F(T)*EEE) >= n_thr ]

con w(T) = peso del espectro de retroceso CEvNS, y (n_F, p_F) los parametros
binomiales del paso NEST (p_F = 1 - F(T), n_F = nint(<N_e>/p_F)).  p_eff = p_F*EEE
exactamente como corre el pipeline hoy (p_vida == 1, no hay vida media modelada).

DOS METODOS INDEPENDIENTES
  (A) ANALITICO -- identidad exacta de la binomial, sin modelo de juguete:
        d/dp Pr[Bin(n,p) >= k] = n * C(n-1, k-1) * p^(k-1) * (1-p)^(n-k)
      y por la regla de la cadena (p = p_F*EEE  =>  dp/dEEE = p_F):
        dS/dEEE|_T = p_F * n_F * C(n_F-1, n_thr-1) * p_eff^(n_thr-1) * (1-p_eff)^(n_F-n_thr)
      w(T), n_F(T), p_F(T) se leen de datos/eee_kappa_grid_<TAG>.dat, que
      chi2_ideal_nest.f90 vuelca con env DUMP_EEE_KAPPA=1 usando EXACTAMENTE la
      misma malla y el mismo paso NEST que el analisis idealizado (Asimov).

  (B) DIFERENCIA FINITA sobre el pipeline compilado:
        kappa ~= ln(S+/S-) / (2*delta) ,  delta = 1%
      - ideal (Xe y Ar): re-corre chi2_ideal con EEE*(1 +- delta) parcheando
        el 'parameter EEE' de constants.f90, y lee R_tot_ROI de
        sensib_ideal_<TAG>.dat por umbral N_e >= n_thr.
      - Xe real (secundario): re-corre red100PE -> chi2 con EEE*(1 +- delta);
        la senal SM aceptada en la ROI de PE (Sum R_SM) recorre la cadena
        ENTERA: binomial(p_F*EEE) + convolucion de plantillas de electron unico
        (mod_detector) + eff_ROI(4:7) + interpolacion a bins de PE.

FUENTES DE LOS VALORES DE EEE (no se inventa ninguno)
  Xe : EEE = 0.328 +- 0.028  (arXiv:2411.18641, l.369 y l.809-812; sigma_rel = 8.54%)
  Ar : EEE = 0.99 ; sin sigma reportada. ReD arXiv:2510.16404 (l.274):
       "The extraction field guarantees 100% extraction efficiency". El paper
       especifico de Ar de RED-100 (ref.[46], Physics 5, 492 (2023)) NO esta en
       el directorio. -> se calcula kappa_Ar, pero el sistematico EEE de Ar se
       declara NO MODELADO (dominante en Ar es Q_y, no EEE).

Salida: stdout (tabla de auditoria) + datos/kappa_eee.dat
SOLO CALCULO. La graficacion va aparte en python/kappa_eee_plot.py.
"""
import os
import re
import shutil
import subprocess
import numpy as np
from scipy.stats import binom
from scipy.special import gammaln

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
DATOS = f"{ROOT}/datos"
DIR = {"Xe": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe",
       "Ar": f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIAr"}
SCRATCH = os.environ.get("KAPPA_SCRATCH", "/tmp/kappa_eee_backup")

EEE_NOM = {"Xe": 0.328, "Ar": 0.99}
EEE_SIGREL = {"Xe": 2.8 / 32.8, "Ar": None}      # Ar: no reportada
NTHR_PHYS = {"Xe": 4, "Ar": 1}                   # piso fisico del ROI de cada blanco
DELTA = 0.01                                     # paso de la diferencia finita

MODS_IDEAL = ("constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 "
              "Tnr_to_e.f90 chi2_ideal_nest.f90")
MODS_PE = ("constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 "
           "Tnr_to_e.f90 mod_detector.f90 red100PE.f90")
MODS_CHI2 = "constants.f90 chi2.f90"
GF = "gfortran -O2 -ffree-line-length-none -o".split()


# --------------------------------------------------------------------------
#  utilidades de compilacion / parcheo de EEE
# --------------------------------------------------------------------------
def build(folder, out, mods):
    subprocess.run(GF + [out] + mods.split(), cwd=folder, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def run(folder, exe):
    return subprocess.run([f"./{exe}"], cwd=folder, check=True,
                          capture_output=True, text=True).stdout


def patch_eee(folder, value):
    """Reemplaza el literal numerico de 'parameter :: EEE = <x>_dp' en constants.f90."""
    path = f"{folder}/constants.f90"
    txt = open(path, encoding="utf-8").read()
    new, n = re.subn(r"(::\s*EEE\s*=\s*)[0-9.eE+]+(_dp)",
                     rf"\g<1>{value:.12f}\g<2>", txt)
    if n != 1:
        raise RuntimeError(f"patch EEE fallo en {path}: {n} coincidencias")
    open(path, "w", encoding="utf-8").write(new)


def read_sensib_rtot(tag):
    """{n_thr: R_tot_ROI}  para F_mode == 'Fnest' de datos/sensib_ideal_<tag>.dat."""
    out = {}
    for ln in open(f"{DATOS}/sensib_ideal_{tag}.dat", encoding="utf-8"):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if p[1] == "Fnest":
            out[int(p[0])] = float(p[2])          # col 2 = R_tot_ROI[ev/kgd]
    return out


def full_chain_S_xe():
    """(Sum R_SM pred ROI [de chi2],  Sum ROI Ne4-7 pre-eff_ROI [de red100PE])."""
    out_pe = run(DIR["Xe"], "red100PE")
    s_pre = float(re.search(r"Sum ROI \(Ne 4-7, antes de eff_ROI\)\s*=\s*([0-9.eE+-]+)",
                            out_pe).group(1))
    out_chi2 = run(DIR["Xe"], "chi2")
    s_fit = float(re.search(r"Sum R_SM \(pred\. ROI\)\s*=\s*([0-9.eE+-]+)",
                            out_chi2).group(1))
    return s_fit, s_pre


# --------------------------------------------------------------------------
#  metodo analitico
# --------------------------------------------------------------------------
def load_grid(tag):
    d = np.loadtxt(f"{DATOS}/eee_kappa_grid_{tag}.dat")
    T, w, nF, pF, lam = d[:, 0], d[:, 1], d[:, 2].astype(int), d[:, 3], d[:, 4]
    return T, w, nF, pF, lam


def log_binom_coeff(n, k):
    """ln C(n,k) vectorizado; -inf donde k<0 o k>n."""
    n = np.asarray(n, float); k = np.asarray(k, float)
    out = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    out[(k < 0) | (k > n)] = -np.inf
    return out


def kappa_analytic(tag, EEE, n_thr):
    """Devuelve (kappa, S, dS, dominante) para el metodo analitico."""
    T, w, nF, pF, lam = load_grid(tag)
    p_eff = pF * EEE
    ok = (nF >= n_thr) & (nF > 0) & (w > 0) & (p_eff > 0) & (p_eff < 1)

    # S_T = Pr[Bin(nF, p_eff) >= n_thr]   (funcion de supervivencia exacta)
    S_T = np.zeros_like(w)
    S_T[ok] = binom.sf(n_thr - 1, nF[ok], p_eff[ok])

    # dS_T/dEEE = p_F * nF * C(nF-1, n_thr-1) * p_eff^(n_thr-1) * (1-p_eff)^(nF-n_thr)
    dS_T = np.zeros_like(w)
    lc = log_binom_coeff(nF[ok] - 1, n_thr - 1)
    log_term = (lc
                + (n_thr - 1) * np.log(p_eff[ok])
                + (nF[ok] - n_thr) * np.log1p(-p_eff[ok]))
    dS_T[ok] = pF[ok] * nF[ok] * np.exp(log_term)

    S = np.sum(w * S_T)
    dS = np.sum(w * dS_T)
    kappa = EEE * dS / S if S > 0 else np.nan

    contrib = w * dS_T
    idom = int(np.argmax(contrib))
    dom = dict(i=idom, T=T[idom], w=w[idom], nF=int(nF[idom]), pF=pF[idom],
               p_eff=p_eff[idom], mean_creados=nF[idom] * pF[idom],
               mean_bin_nm1=(nF[idom] - 1) * p_eff[idom],
               frac_contrib=contrib[idom] / np.sum(contrib))
    return kappa, S, dS, dom, (T, w, nF, pF, p_eff, S_T, dS_T)


# --------------------------------------------------------------------------
#  metodo de diferencia finita
# --------------------------------------------------------------------------
def kappa_fd_ideal(tag):
    """{n_thr: kappa_fd} re-corriendo chi2_ideal con EEE*(1 +- delta)."""
    folder = DIR[tag]
    Sp = Sm = None
    try:
        patch_eee(folder, EEE_NOM[tag] * (1 + DELTA))
        build(folder, "chi2_ideal", MODS_IDEAL)
        run(folder, "chi2_ideal")
        Sp = read_sensib_rtot(tag)
        patch_eee(folder, EEE_NOM[tag] * (1 - DELTA))
        build(folder, "chi2_ideal", MODS_IDEAL)
        run(folder, "chi2_ideal")
        Sm = read_sensib_rtot(tag)
    finally:
        shutil.copy(f"{SCRATCH}/constants_{tag}.f90", f"{folder}/constants.f90")
        build(folder, "chi2_ideal", MODS_IDEAL)
        run(folder, "chi2_ideal")                     # deja datos/ en nominal
    return {k: np.log(Sp[k] / Sm[k]) / (2 * DELTA) for k in Sp}


def kappa_fd_fullchain_xe():
    """kappa_fd sobre red100PE -> chi2 (cadena completa con SE + eff_ROI)."""
    folder = DIR["Xe"]
    try:
        patch_eee(folder, EEE_NOM["Xe"] * (1 + DELTA))
        build(folder, "red100PE", MODS_PE)
        build(folder, "chi2", MODS_CHI2)
        Sp_fit, Sp_pre = full_chain_S_xe()
        patch_eee(folder, EEE_NOM["Xe"] * (1 - DELTA))
        build(folder, "red100PE", MODS_PE)
        build(folder, "chi2", MODS_CHI2)
        Sm_fit, Sm_pre = full_chain_S_xe()
    finally:
        shutil.copy(f"{SCRATCH}/constants_Xe.f90", f"{folder}/constants.f90")
        build(folder, "red100PE", MODS_PE)
        build(folder, "chi2", MODS_CHI2)
        full_chain_S_xe()
    return (np.log(Sp_fit / Sm_fit) / (2 * DELTA),
            np.log(Sp_pre / Sm_pre) / (2 * DELTA))


# --------------------------------------------------------------------------
#  main
# --------------------------------------------------------------------------
def main():
    os.makedirs(SCRATCH, exist_ok=True)
    for tag in ("Xe", "Ar"):
        shutil.copy(f"{DIR[tag]}/constants.f90", f"{SCRATCH}/constants_{tag}.f90")

    # regenerar los grids desde la fuente actual (por si acaso), sin tocar nada mas
    for tag in ("Xe", "Ar"):
        build(DIR[tag], "chi2_ideal", MODS_IDEAL)
        subprocess.run(["./chi2_ideal"], cwd=DIR[tag],
                       env={**os.environ, "DUMP_EEE_KAPPA": "1"},
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    lines = []
    def emit(s=""):
        print(s); lines.append(s)

    emit("=" * 92)
    emit(" kappa_EEE = (EEE/S) dS/dEEE   --   amplificacion de la incertidumbre de EEE")
    emit("=" * 92)
    emit(f" EEE_nom:  Xe = {EEE_NOM['Xe']}  (sigma_rel = {EEE_SIGREL['Xe']*100:.2f}%,"
         f" arXiv:2411.18641)")
    emit(f"           Ar = {EEE_NOM['Ar']}  (sigma NO reportada; ReD dice extraccion"
         f" ~100% por diseno de campo)")
    emit(f" delta (diferencia finita) = {DELTA*100:.1f}%   |   malla T: n_T=800"
         f" (identica a chi2_ideal_nest)")

    # ---- analitico: barrido de umbral + auditoria del termino dominante ----
    ana = {}
    for tag in ("Xe", "Ar"):
        emit()
        emit("-" * 92)
        emit(f" [{tag}]  METODO ANALITICO  (identidad exacta de la binomial)")
        emit("-" * 92)
        T, w, nF, pF, lam = load_grid(tag)
        wmean_T = np.sum(w * T) / np.sum(w)
        emit(f"  nodos de malla con senal (n_F>0): {int(np.sum(nF > 0))}/800   "
             f"| <T>_w = {wmean_T:.3f} keV   | p_F rango = "
             f"[{pF[nF > 0].min():.3f}, {pF[nF > 0].max():.3f}]")

        emit(f"  P0(n_e0) = Binomial(n_F(T), p_F(T)) -- primeros 10 n_e0, en T dominante:")
        ana[tag] = {}
        for nthr in (1, 2, 3, 4):
            k, S, dS, dom, _ = kappa_analytic(tag, EEE_NOM[tag], nthr)
            ana[tag][nthr] = (k, S, dS, dom)

        dom = ana[tag][NTHR_PHYS[tag]][3]
        p0 = binom.pmf(np.arange(10), dom["nF"], dom["pF"])
        emit(f"    T_dom = {dom['T']:.4f} keV   n_F = {dom['nF']}   p_F = {dom['pF']:.4f}"
             f"   p_eff = p_F*EEE = {dom['p_eff']:.4f}")
        emit("    n_e0 :  " + "  ".join(f"{i:8d}" for i in range(10)))
        emit("    P0   :  " + "  ".join(f"{v:8.5f}" for v in p0))
        emit(f"    <N_e creados> = n_F*p_F = {dom['mean_creados']:.3f} ;  "
             f"media de Bin(n_F-1, p_eff) = {dom['mean_bin_nm1']:.3f}")

        emit()
        emit(f"  {'n_thr':>6} {'S (int w*S_T)':>16} {'dS/dEEE':>16} {'kappa_EEE':>12}"
             f"   {'T_dom[keV]':>10} {'n_F_dom':>8} {'p_eff_dom':>10}"
             f" {'frac_dom':>9}   interpretacion")
        for nthr in (1, 2, 3, 4):
            k, S, dS, dom = ana[tag][nthr]
            note = ""
            if nthr == NTHR_PHYS[tag]:
                note = "<-- piso fisico del ROI"
            emit(f"  {nthr:>6} {S:>16.6e} {dS:>16.6e} {k:>12.4f}"
                 f"   {dom['T']:>10.4f} {dom['nF']:>8d} {dom['p_eff']:>10.4f}"
                 f" {dom['frac_contrib']:>9.3f}   {note}")

    # ---- diferencia finita: ideal (Xe y Ar) ----
    emit()
    emit("-" * 92)
    emit(" METODO DE DIFERENCIA FINITA  (pipeline compilado, delta = 1%)")
    emit("-" * 92)
    fd_ideal = {}
    for tag in ("Xe", "Ar"):
        fd_ideal[tag] = kappa_fd_ideal(tag)
        emit(f"  [{tag}] ideal  (S = R_tot_ROI de sensib_ideal_{tag}.dat, F de NEST):")
        for nthr in sorted(fd_ideal[tag]):
            ka = ana[tag][nthr][0]
            kf = fd_ideal[tag][nthr]
            rel = abs(kf - ka) / abs(ka) if ka else float("nan")
            emit(f"     n_thr={nthr}:  kappa_analitico = {ka:>9.4f}   "
                 f"kappa_dif.finita = {kf:>9.4f}   |dif| rel = {rel:6.2%}")

    # ---- diferencia finita: Xe cadena completa ----
    kf_fit, kf_pre = kappa_fd_fullchain_xe()
    emit()
    emit("  [Xe] cadena COMPLETA red100PE -> chi2  (binomial + convolucion SE + eff_ROI):")
    emit(f"     S = 'Sum R_SM (pred. ROI)' [entra al ajuste chi2] :  kappa = {kf_fit:.4f}")
    emit(f"     S = 'Sum ROI Ne4-7' [pre eff_ROI, control]        :  kappa = {kf_pre:.4f}")
    emit(f"     (comparar con ideal N_e>=4:  analitico "
         f"{ana['Xe'][4][0]:.4f} / dif.finita {fd_ideal['Xe'][4]:.4f})")

    # ---- resumen: efecto sobre la senal ----
    emit()
    emit("=" * 92)
    emit(" RESULTADO  (kappa * sigma_rel(EEE) = incertidumbre relativa inducida en la senal)")
    emit("=" * 92)
    for tag in ("Xe", "Ar"):
        nthr = NTHR_PHYS[tag]
        ka = ana[tag][nthr][0]
        kf = fd_ideal[tag][nthr]
        sr = EEE_SIGREL[tag]
        emit(f"  [{tag}] n_thr={nthr} (ideal):  kappa = {ka:.3f} (analitico) /"
             f" {kf:.3f} (dif.finita)")
        if sr is None:
            emit(f"        sigma_rel(EEE) NO reportada para Ar -> sistematico EEE "
                 f"NO modelado (dominante en Ar: Q_y).")
        else:
            emit(f"        sigma_rel(EEE) = {sr:.2%}  ->  "
                 f"kappa*sigma_rel = {ka*sr:.2%} de incertidumbre relativa en S.")
    emit()
    emit(f"  [Xe] cadena completa (real, N_e 4-7 + SE + eff_ROI):  kappa = {kf_fit:.3f}"
         f"  ->  kappa*sigma_rel = {kf_fit*EEE_SIGREL['Xe']:.2%}")

    with open(f"{DATOS}/kappa_eee.dat", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n  -> {DATOS}/kappa_eee.dat")


if __name__ == "__main__":
    main()
