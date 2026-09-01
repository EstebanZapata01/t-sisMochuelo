#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VALIDACION del motor NSI contra CONUS+ (De Romeri, Papoulias, Sanchez Garcia,
Phys. Rev. D 111, 075025 (2025)).

A diferencia de RED-100 (Xe: limite superior debil; Ar: proyeccion sin datos),
CONUS+ TIENE datos reales con exceso a 3.7 sigma y publica cotas numericas
(Tabla II). Reproducirlas con el mismo motor chi2 (analitico en alpha) que se
uso para Xe es una validacion independiente de TODA la maquinaria NSI:
prediccion (resolucion Ge + quenching Lindhard + seccion eficaz) + estadistica.

Metodo:
  1. Se corre FORTRAN90/N_EventosCEvNS_NSI/chi2_nsi (2pchi2.f90) que escribe
     datos/chi2_nsi_2Dconus.dat  (grilla eps_x eps_y chi2, ipar fijo).
  2. Aqui se toma el corte 1D a lo largo de cada eje (el otro parametro = 0),
     se calcula Dchi2 = chi2 - min(chi2) y se extraen los intervalos con
     Dchi2 <= 1  (1 sigma, 1 g.d.l., que es como esta hecha la Tabla II).
  3. Se comparan con los valores publicados.

Uso:
    python3 valida_conus.py                 # usa el .dat ya generado
    (para otro ipar: editar 'ipar' en 2pchi2.f90, recompilar y re-correr)
"""
import numpy as np
import os

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"

# --- Tabla II de Phys. Rev. D 111, 075025 (CONUS+, "This Work", 1 sigma) ------
#     (analisis CEvNS-only; cada parametro con los demas puestos a cero)
PUBLICADO = {
    r"$\epsilon_{ee}^{dV}$":   [(-0.034, 0.024), (0.322, 0.380)],
    r"$\epsilon_{ee}^{uV}$":   [(-0.037, 0.026), (0.348, 0.411)],
    r"$\epsilon_{e\mu}^{dV}$": [(-0.114, 0.114)],
    r"$\epsilon_{e\mu}^{uV}$": [(-0.123, 0.123)],
    r"$\epsilon_{e\tau}^{dV}$": [(-0.114, 0.114)],
    r"$\epsilon_{e\tau}^{uV}$": [(-0.123, 0.123)],
}


def carga_grilla(path):
    d = np.loadtxt(path)
    n = int(round(np.sqrt(d.shape[0])))
    X = d[:, 0].reshape(n, n)
    Y = d[:, 1].reshape(n, n)
    Z = d[:, 2].reshape(n, n)
    return X, Y, Z


def etiquetas(path):
    with open(path) as f:
        ln = [l.strip() for l in f if l.strip()]
    return ln[0], ln[1]


def intervalos_por_debajo(x, dchi2, umbral):
    """Devuelve la lista de (x_lo, x_hi) donde dchi2(x) <= umbral,
    con interpolacion lineal en los cruces."""
    dentro = dchi2 <= umbral
    tramos = []
    i = 0
    n = len(x)
    while i < n:
        if not dentro[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and dentro[j + 1]:
            j += 1
        # borde inferior
        if i == 0:
            xlo = x[0]
        else:
            f = (umbral - dchi2[i - 1]) / (dchi2[i] - dchi2[i - 1])
            xlo = x[i - 1] + f * (x[i] - x[i - 1])
        # borde superior
        if j == n - 1:
            xhi = x[-1]
        else:
            f = (umbral - dchi2[j]) / (dchi2[j + 1] - dchi2[j])
            xhi = x[j] + f * (x[j + 1] - x[j])
        tramos.append((xlo, xhi))
        i = j + 1
    return tramos


def corte_1d(X, Y, Z, eje):
    """Corte a lo largo del eje pedido con el otro parametro = 0.
       eje = 'x'  -> varia X, fija Y=0 (fila mas cercana a Y=0)
       eje = 'y'  -> varia Y, fija X=0 (columna mas cercana a X=0)"""
    if eje == "x":
        k = int(np.argmin(np.abs(Y[:, 0])))
        return X[k, :].copy(), Z[k, :].copy()
    else:
        k = int(np.argmin(np.abs(X[0, :])))
        return Y[:, k].copy(), Z[:, k].copy()


def fmt_int(tramos):
    return "  U  ".join(f"[{a:+.3f}, {b:+.3f}]" for a, b in tramos)


def compara(nombre, tramos_sim, umbral):
    pub = PUBLICADO.get(nombre)
    print(f"\n  {nombre}   (Dchi2 <= {umbral})")
    print(f"    simulacion : {fmt_int(tramos_sim)}")
    if pub is None:
        print("    publicado  : (no en la tabla)")
        return
    print(f"    publicado  : {fmt_int(pub)}")
    # comparacion cruda: casar tramo a tramo por cercania de centro
    if len(tramos_sim) == len(pub):
        difs = []
        for (a, b), (pa, pb) in zip(sorted(tramos_sim), sorted(pub)):
            difs.append(abs(a - pa))
            difs.append(abs(b - pb))
        print(f"    |dif| bordes: max {max(difs):.3f} , media {np.mean(difs):.3f}")
    else:
        print(f"    (num. de tramos distinto: sim {len(tramos_sim)} vs pub {len(pub)})")


def main():
    f_dat = f"{BASE}/chi2_nsi_2Dconus.dat"
    f_cfg = f"{BASE}/nsi_config_conus.txt"
    if not os.path.exists(f_dat):
        raise SystemExit(f"Falta {f_dat}. Compila y corre FORTRAN90/"
                         "N_EventosCEvNS_NSI/chi2_nsi primero.")
    X, Y, Z = carga_grilla(f_dat)
    xl, yl = etiquetas(f_cfg)

    print("=" * 74)
    print(" VALIDACION del motor NSI contra CONUS+  (Phys. Rev. D 111, 075025)")
    print(f" grilla {Z.shape[0]}x{Z.shape[1]} , plano leido de nsi_config_conus.txt:")
    print(f"   eje x = {xl}    eje y = {yl}")
    print("=" * 74)

    UMBRAL_1SIG = 1.0        # 1 g.d.l., como la Tabla II
    UMBRAL_90 = 2.706        # 1 g.d.l., 90 % C.L. (referencia)

    for eje, lab in (("x", xl), ("y", yl)):
        x, z = corte_1d(X, Y, Z, eje)
        # ordenar por si el eje viene invertido
        o = np.argsort(x)
        x, z = x[o], z[o]
        dchi2 = z - z.min()
        xbest = x[np.argmin(z)]
        print(f"\n  --- corte 1D en {lab}  (otro parametro = 0) ---")
        print(f"      chi2_min en {lab} = {xbest:+.3f}   (chi2_min = {z.min():.3f})")
        compara(lab, intervalos_por_debajo(x, dchi2, UMBRAL_1SIG), UMBRAL_1SIG)
        tr90 = intervalos_por_debajo(x, dchi2, UMBRAL_90)
        print(f"      (90 % C.L., Dchi2<=2.706): {fmt_int(tr90)}")

    # region 2D permitida al 90 % (2 g.d.l.) -> comparacion cualitativa con Fig. 4
    d2 = Z - Z.min()
    frac90 = float((d2 <= 4.605).mean())
    print(f"\n  --- region 2D ---")
    print(f"      chi2_min global = {Z.min():.3f}")
    print(f"      fraccion de la caja |eps|<=1 con Dchi2<=4.605 (90%, 2 gdl) = {frac90:.4f}")
    print(f"      (Fig. 4 del paper: dos bandas / un anillo, NO toda la caja)")
    print("=" * 74)


if __name__ == "__main__":
    main()
