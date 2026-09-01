#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparacion IDEAL SIMETRICA Xe vs Ar de la SENSIBILIDAD PROYECTADA de RED-100
al CEvNS de antineutrinos de reactor y a NSI.

*** Todo aqui es una SENSIBILIDAD (proyeccion Asimov), NO una exclusion con
    datos. Se asume: sustraccion perfecta del fondo, cero sistematicos,
    eficiencia de seleccion = 1. Es el LIMITE OPTIMISTA de cada blanco.
    El unico numero con datos reales es el del Xe (rama chi2.f90,
    A_90 ~ 111 xSM), que se reporta aparte como referencia. ***

Entradas (salidas de FORTRAN90/N_EventosCEvNS_NSI{Xe,Ar}/chi2_ideal_nest.f90):
    datos/espectro_Ne_ideal_{Xe,Ar}.dat
    datos/sensib_ideal_{Xe,Ar}.dat
    datos/chi2_nsi_2D{Xe,Ar}_ideal.dat
    datos/nsi_config_ideal_{Xe,Ar}.txt   (etiquetas de ejes segun `ipar`)

Salidas (datos/), en escala de grises, serif, una figura por concepto:
    fig_ideal_1_espectro_Ne.png
    fig_ideal_2_retencion_umbral.png
    fig_ideal_3_sensibilidad_exposicion.png
    fig_ideal_4_figura_merito.png
    fig_ideal_5_plano_NSI.png
    fig_ideal_XeAr.png            (panel resumen 2x2, opcional)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Circle

BASE = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos"
MULT = np.array([1, 2, 5, 10, 50, 100, 335], dtype=float)
EXPO_BASE = 192.0
XE_REAL_A90 = 111.1
COH_ONLY = 3.8

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "axes.linewidth": 0.9, "axes.grid": True,
    "grid.color": "0.82", "grid.linewidth": 0.5,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.fontsize": 9,
    "lines.linewidth": 1.7,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 220, "savefig.bbox": "tight",
})

# --------------------------------------------------------------- lectura
def load_sensib(tag):
    out = {}
    for ln in open(f"{BASE}/sensib_ideal_{tag}.dat"):
        if ln.lstrip().startswith("#") or not ln.strip():
            continue
        p = ln.split()
        out[(int(p[0]), p[1])] = dict(
            R_tot=float(p[2]), frac=float(p[3]), Ntot192=float(p[4]),
            A90=np.array([float(x) for x in p[5:12]]))
    return out

def load_spec(tag):
    d = np.loadtxt(f"{BASE}/espectro_Ne_ideal_{tag}.dat")
    return d[:, 0], d[:, 1], d[:, 2]

def load_nsi(tag):
    d = np.loadtxt(f"{BASE}/chi2_nsi_2D{tag}_ideal.dat")
    n = int(round(np.sqrt(d.shape[0])))
    return d[:, 0].reshape(n, n), d[:, 1].reshape(n, n), d[:, 2].reshape(n, n)

def load_nsi_labels(tag):
    """Etiqueta LaTeX de los ejes del plano NSI, escrita por chi2_ideal_nest.f90
    segun `ipar`. Si falta el archivo, cae al plano por defecto (ipar=5)."""
    fallback = (r"$\varepsilon_{ee}^{dV}$", r"$\varepsilon_{e\mu}^{dV}$")
    try:
        raw = [ln.strip() for ln in open(f"{BASE}/nsi_config_ideal_{tag}.txt")
               if ln.strip()]
        xl, yl = raw[0], raw[1]
    except (FileNotFoundError, IndexError):
        return fallback
    # el codigo Fortran escribe \epsilon; el resto del script usa \varepsilon
    return xl.replace(r"\epsilon", r"\varepsilon"), yl.replace(r"\epsilon", r"\varepsilon")

def blind_spot_from_grid(X, Y, Z):
    """Punto ciego (Dchi2 -> 0) sobre el eje eps_y = 0, leido de la grilla:
    no se hardcodea, vale para cualquier `ipar` que tenga un cero en ese eje."""
    jrow = int(np.argmin(np.abs(Y[:, 0])))          # fila mas cercana a eps_y = 0
    zr = Z[jrow]
    xr = X[jrow]
    cand = xr[(zr < 1.0) & (xr > 1e-3)]             # cerca de Dchi2=0, lado eps>0
    if cand.size:
        return float(xr[np.argmin(np.where(xr > 1e-3, zr, np.inf))])
    return float("nan")

sx, sa = load_sensib("Xe"), load_sensib("Ar")
ne, rx, _ = load_spec("Xe")
_,  ra, _ = load_spec("Ar")
THR = np.array([1, 2, 3, 4])
expo = EXPO_BASE * MULT

# --------------------------------------------------------------- helper figura+caption
def new_fig(w=7.8, h=8.0, ratio=(2.25, 1.25)):
    fig, (ax, cx) = plt.subplots(
        2, 1, figsize=(w, h), layout="constrained",
        gridspec_kw=dict(height_ratios=list(ratio)))
    fig.get_layout_engine().set(h_pad=0.14, hspace=0.10)
    cx.axis("off")
    cx.set_xlim(0, 1); cx.set_ylim(0, 1)
    return fig, ax, cx

def put_caption(cx, head, body):
    cx.add_patch(plt.Rectangle((0.005, 0.02), 0.99, 0.96, transform=cx.transAxes,
                               fill=False, ec="0.75", lw=0.7))
    cx.text(0.03, 0.90, head, transform=cx.transAxes, va="top", ha="left",
            fontsize=9, fontweight="bold")
    cx.text(0.03, 0.72, body, transform=cx.transAxes, va="top", ha="left",
            fontsize=8, linespacing=1.55)

def save(fig, name):
    out = f"{BASE}/{name}"
    fig.savefig(out)
    plt.close(fig)
    print(f"  {out}")

# estilos Xe / Ar consistentes
def xe_plot(ax, x, y, **kw):
    kw = {**dict(color="black", ls="-", marker="o", mfc="white", mec="black", ms=5.5), **kw}
    return ax.plot(x, y, **kw)

def ar_plot(ax, x, y, **kw):
    kw = {**dict(color="0.45", ls="--", marker="s", mfc="0.45", mec="0.45", ms=5), **kw}
    return ax.plot(x, y, **kw)

print("Generando figuras separadas...")

# =====================================================================
# 1) ESPECTRO DE IONIZACION
# =====================================================================
fig, ax, cx = new_fig()
ax.step(ne, rx, where="mid", color="black", lw=1.7)
ax.step(ne, ra, where="mid", color="0.45", ls="--", lw=1.7)
xe_plot(ax, ne, rx, ls="none", label="Xe (LXe)")
ar_plot(ax, ne, ra, ls="none", label="Ar (LAr)")
for t in (2, 3, 4):
    ax.axvline(t - 0.5, color="0.6", lw=0.8, ls=":")
ax.set_yscale("log"); ax.set_xlim(0.4, 21); ax.set_ylim(3e-6, 40)
ax.set_xlabel(r"$N_e$   (electrones de ionización extraídos)")
ax.set_ylabel(r"tasa CEvNS SM   $R(N_e)$   [ev / (kg$\cdot$día)]")
ax.set_title("Espectro de ionización  —  Xe y Ar por el MISMO código")
ax.legend(loc="upper right")
ax.annotate("Xe: ~80 % de la señal\nen $N_e\\,{=}\\,1$; luego se hunde",
            xy=(1.0, rx[0]), xytext=(4.6, 3.2), fontsize=8.2, color="0.15",
            arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8,
                            connectionstyle="arc3,rad=-0.25"))
ax.annotate("Ar: decae suave,\ncola hasta $N_e\\sim20$",
            xy=(13, ra[12]), xytext=(9.5, 6e-3), fontsize=8.2, color="0.15",
            arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8))
put_caption(cx, "Qué muestra",
    "Tasa de eventos CEvNS del Modelo Estándar por número de electrones de ionización,\n"
    "con el MISMO programa (chi2_ideal_nest.f90) para los dos blancos: sólo cambian $Q_W$,\n"
    "la masa nuclear y la tabla de yield de NEST.  Punteadas: umbrales $N_e\\!\\geq\\!2,3,4$.\n"
    "\n"
    "Lectura: los retrocesos de Xe son blandos ($T_{max}\\!\\sim\\!1$ keV) y su yield es $\\sim$0 en\n"
    "el umbral $\\Rightarrow$ ~80 % de la señal cae en $N_e\\!=\\!1$.  Los de Ar son más duros\n"
    "($T_{max}\\!\\sim\\!3.4$ keV) $\\Rightarrow$ el espectro se reparte hasta $N_e\\!\\sim\\!20$.  Esta forma\n"
    "explica todos los resultados siguientes.")
save(fig, "fig_ideal_1_espectro_Ne.png")

# =====================================================================
# 2) RETENCION DE SENAL vs UMBRAL
# =====================================================================
fig, ax, cx = new_fig()
fxv = np.array([sx[(t, "Fnest")]["frac"] for t in THR]) * 100
fav = np.array([sa[(t, "Fnest")]["frac"] for t in THR]) * 100
xe_plot(ax, THR, fxv, ms=7); ar_plot(ax, THR, fav, ms=6.5)
for t, v in zip(THR, fxv):
    ax.annotate(f"{v:.1f} %", (t, v), textcoords="offset points", xytext=(0, 10),
                ha="center", fontsize=8.5)
for t, v in zip(THR, fav):
    ax.annotate(f"{v:.0f} %", (t, v), textcoords="offset points", xytext=(0, -15),
                ha="center", fontsize=8.5, color="0.30")
ax.set_yscale("log"); ax.set_xticks(THR); ax.set_xlim(0.75, 4.25); ax.set_ylim(0.3, 260)
ax.set_xlabel(r"umbral inferior aplicado   $N_e \geq$")
ax.set_ylabel(r"señal CEvNS retenida   [% de la de $N_e\!\geq\!1$]")
ax.set_title("Coste de subir el umbral inferior")
ax.legend(["Xe (LXe)", "Ar (LAr)"], loc="lower left")
put_caption(cx, "Qué muestra",
    "Fracción de la tasa CEvNS total que sobrevive al exigir un umbral inferior en $N_e$,\n"
    "relativa a la ventana más amplia $N_e\\geq1$.  El umbral $N_e\\geq4$ es el de RED-100-Xe,\n"
    "impuesto por el fondo de electrón único (ruido instrumental), NO por sensibilidad.\n"
    "Lectura: subir a $N_e\\geq4$ deja al Xe con el 0.6 % de su señal (cae en $N_e\\!=\\!1$),\n"
    "mientras el Ar conserva el 40 %.  Por eso comparar Xe@4..7 con Ar@1..5 daba un\n"
    "factor artificial ~470 a favor de Ar: era la ventana, no la física.")
save(fig, "fig_ideal_2_retencion_umbral.png")

# =====================================================================
# 3) SENSIBILIDAD PROYECTADA vs EXPOSICION
# =====================================================================
fig, ax, cx = new_fig()
xe_plot(ax, expo, sx[(1, "Fnest")]["A90"], label=r"Xe  proyección ($N_e\!\geq\!1$)")
ar_plot(ax, expo, sa[(1, "Fnest")]["A90"], label=r"Ar  proyección ($N_e\!\geq\!1$)")
ax.plot(expo, sx[(4, "Fnest")]["A90"], color="black", ls=":", marker="^",
        mfc="black", mec="black", ms=5, lw=1.1, label=r"Xe  proyección ($N_e\!\geq\!4$)")
ax.axhline(XE_REAL_A90, color="black", lw=1.3, ls=(0, (7, 3)))
ax.text(expo[0] * 1.08, XE_REAL_A90 * 0.60,
        f"Xe con DATOS reales (ajuste ON$-$OFF 2024)  $\\approx$ {XE_REAL_A90:.0f} $\\times$ SM",
        fontsize=8.6, va="top")
ax.annotate("", xy=(expo[3], XE_REAL_A90 * 0.78), xytext=(expo[3], 1.7),
            arrowprops=dict(arrowstyle="<->", color="0.35", lw=1.1))
ax.text(expo[3] * 1.18, 13,
        "brecha $\\approx\\times100$:\nfondo real + cortes\n+ sistemáticos\n+ ajuste a 1 histograma",
        fontsize=8.1, color="0.15", va="center")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.95, 270)
ax.set_xlabel(r"exposición reactor ON   [kg$\cdot$día]      (192 = dato 2024;  64 000 $\approx$ 1 año)")
ax.set_ylabel(r"$A_{90}$   [amplitud CEvNS  $\times$ SM,  90 % C.L.]")
ax.set_title("Sensibilidad PROYECTADA a la amplitud CEvNS")
ax.legend(loc="center left")
ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
ax.yaxis.set_minor_formatter(NullFormatter())
put_caption(cx, "Qué muestra — y por qué es un LÍMITE DE SENSIBILIDAD, no una medida",
    "$A_{90}$ = cuánto tendría que desviarse la amplitud CEvNS ($\\propto Q_W^2$) del SM para\n"
    "ser detectada al 90 % C.L.  Las curvas son una proyección Asimov de estadística\n"
    "pura (sin fondo, sin sistemáticos, eff = 1): el MEJOR caso posible, no lo alcanzable.\n"
    "La línea a trazos es el único número con datos reales (RED-100 con Xe).\n"
    "\n"
    "Lectura: la proyección ideal roza $A_{90}\\!\\approx\\!1.02$–1.04; el Xe real está en 111,\n"
    "$\\sim\\!100\\times$ peor.  Por eso la comparación Xe vs Ar sólo es legítima entre las dos\n"
    "proyecciones ideales, NUNCA contra el 111.")
save(fig, "fig_ideal_3_sensibilidad_exposicion.png")

# =====================================================================
# 4) FIGURA DE MERITO INTRINSECA
# =====================================================================
fig, ax, cx = new_fig()
rn = np.array([sa[(t, "Fnest")]["R_tot"] / sx[(t, "Fnest")]["R_tot"] for t in THR])
r1 = np.array([sa[(t, "F1")]["R_tot"] / sx[(t, "F1")]["R_tot"] for t in THR])
ax.plot(THR, rn, color="black", ls="-", marker="D", mfc="white", mec="black",
        ms=6.5, label="F de NEST (cada blanco)")
ax.plot(THR, r1, color="0.5", ls="--", marker="D", mfc="0.5", mec="0.5",
        ms=5.5, label=r"$F=1$ (Poisson) en ambos")
ax.axhline(1.0, color="0.6", lw=0.8)
ax.axhline(COH_ONLY, color="black", lw=1.0, ls=":")
ax.text(2.5, COH_ONLY * 1.13,
        r"predicción de la coherencia sola:  Xe $\times3.8$  (Ar/Xe = 1/3.8 $\approx$ 0.26)",
        fontsize=7.9, ha="center")
for t, v in zip(THR, rn):
    ax.annotate(f"{v:.0f}" if v >= 10 else f"{v:.1f}", (t, v),
                textcoords="offset points", xytext=(8, -2), fontsize=8.5, fontweight="bold")
ax.set_yscale("log"); ax.set_xticks(THR); ax.set_xlim(0.75, 4.4); ax.set_ylim(0.7, 400)
ax.set_xlabel(r"umbral inferior   $N_e \geq$")
ax.set_ylabel(r"eventos CEvNS detectados   Ar / Xe")
ax.set_title("Figura de mérito intrínseca:  ¿cuántos eventos ve cada blanco?")
ax.legend(loc="upper left")
put_caption(cx, "Qué muestra",
    "Cociente del número total de eventos CEvNS detectados (Ar entre Xe), con la misma\n"
    "exposición y ventana, en función del umbral.  Como $A_{90}\\!-\\!1 \\propto 1/\\sqrt{N}$, este\n"
    "cociente ES la ventaja relativa en sensibilidad.\n"
    "Lectura: la coherencia sola ($Q_W^2\\times$átomos/kg) predice Xe $\\times3.8$.  Se INVIERTE al\n"
    "exigir señal detectable: en $N_e\\geq1$ el Ar ya ve $\\times3.4$ más eventos, y la ventaja\n"
    "crece con el umbral (hasta $\\times243$ en $N_e\\geq4$) porque el Xe se queda sin espectro.\n"
    "F de NEST vs F=1: el cociente casi no cambia $\\Rightarrow$ el resultado no depende del modelo de fluctuación.")
save(fig, "fig_ideal_4_figura_merito.png")

# =====================================================================
# 5) PLANO NSI 2D  (sensibilidad proyectada)
# =====================================================================
fig, axs = plt.subplots(3, 2, figsize=(9.8, 8.2), layout="constrained",
                        gridspec_kw=dict(height_ratios=[2.35, 0.32, 1.15]))
fig.get_layout_engine().set(h_pad=0.10, hspace=0.06)
axX, axA = axs[0, 0], axs[0, 1]
gs = axs[0, 0].get_gridspec()
for r in (1, 2):
    axs[r, 0].remove(); axs[r, 1].remove()
lx = fig.add_subplot(gs[1, :]); lx.axis("off")
cx = fig.add_subplot(gs[2, :])
cx.axis("off"); cx.set_xlim(0, 1); cx.set_ylim(0, 1)
ZOOM_X, ZOOM_Y = (-0.15, 0.55), (-0.35, 0.35)
info = {}
blind_info = {}
for axi, tag in ((axX, "Xe"), (axA, "Ar")):
    X, Y, Z = load_nsi(tag)
    xl, yl = load_nsi_labels(tag)
    blind = blind_spot_from_grid(X, Y, Z)
    blind_info[tag] = blind
    excl = 100.0 * (Z >= 4.605).mean()
    info[tag] = excl
    axi.contourf(X, Y, Z, levels=[4.605, Z.max() + 1], colors=["0.86"])
    axi.contourf(X, Y, Z, levels=[0.0, 4.605], colors=["white"])
    axi.contour(X, Y, Z, levels=[4.605], colors="black", linewidths=1.3)
    axi.plot(0, 0, "+", color="black", ms=12, mew=1.8)
    axi.annotate("SM", (0, 0), textcoords="offset points", xytext=(-5, 8),
                 fontsize=9, ha="right")
    if np.isfinite(blind):
        # locus ciego q_eff^2 = Q_W^2: en 2D es un circulo que pasa por el SM y
        # corta el eje eps_y=0 en `blind`; centro (blind/2, 0), radio blind/2.
        rc = blind / 2.0
        axi.add_patch(Circle((rc, 0.0), rc, fill=False, ls=(0, (4, 3)),
                             ec="0.35", lw=1.1))
        axi.plot(blind, 0, "x", color="black", ms=8, mew=1.7)
        axi.annotate(f"{xl[:-1]}\\,{{=}}\\,{blind:.2f}$", (blind, 0),
                     textcoords="offset points", xytext=(7, -14), fontsize=8.2)
    axi.set_title(f"{tag}  —  proyección sensible al {excl:.1f} % de $|\\varepsilon|\\!\\leq\\!1$")
    axi.set_xlabel(xl)
    axi.set_ylabel(yl)
    axi.set_xlim(*ZOOM_X); axi.set_ylim(*ZOOM_Y); axi.set_aspect("equal")

handles = [
    Patch(fc="white", ec="black", lw=1.3, label=r"dentro del alcance ($\Delta\chi^2\!<\!4.605$)"),
    Patch(fc="0.86", ec="none", label=r"fuera del alcance (se distinguiría del SM, 90 % C.L.)"),
    Line2D([0], [0], color="black", marker="+", ls="none", mew=1.8, label=r"SM ($\varepsilon=0$)"),
    Line2D([0], [0], color="0.35", ls=(0, (4, 3)), lw=1.1,
           label=r"locus ciego $q_{\rm eff}^2\!=\!Q_W^2$  ($\times$ = corte con $\varepsilon_y\!=\!0$)"),
]
lx.legend(handles=handles, loc="center", ncol=2, fontsize=8.5,
          handletextpad=0.6, columnspacing=1.6, borderaxespad=0.0)
put_caption(cx,
    "Por qué el alcance cubre casi todo el plano — y por qué es SÓLO una proyección",
    "El CEvNS mide la tasa total; la NSI de quark down entra multiplicada por $Z\\!+\\!2N$\n"
    "($\\approx$ 208 en Xe, 62 en Ar) y compite con $Q_W\\!\\approx\\!-37$ / $-11$.  Con ese brazo de\n"
    "palanca un $|\\varepsilon|\\!\\sim\\!0.05$ ya cambia la tasa $>2\\times$: con miles de eventos Asimov y\n"
    "sin fondo, casi cualquier $\\varepsilon$ se separa del SM.  Sólo sobrevive el anillo 'punto\n"
    "ciego', donde $q_{\\rm eff}^2$ vuelve a valer $Q_W^2$ y la tasa es idéntica al SM.\n"
    "CLAVE: NO es una exclusión, es sensibilidad de estadística pura (fondo perfecto,\n"
    "cero sistemáticos).  El análisis REAL de Xe ($A_{90}\\!\\approx\\!111$) no separa del SM\n"
    "ningún punto de esta caja: el alcance real está mucho más cerca de 'nada'.")
save(fig, "fig_ideal_5_plano_NSI.png")

# =====================================================================
# tabla por consola  +  chequeo analitico del alcance NSI
# =====================================================================
L = "-" * 80
print("\n" + "=" * 80)
print(" SENSIBILIDAD PROYECTADA IDEAL  Xe vs Ar   (Asimov conteo puro, 192 kg*dia)")
print("=" * 80)
print(f"{'F':>6} {'N_e>=':>6} | {'R_tot Xe':>10} {'R_tot Ar':>10} | {'Ar/Xe':>7} "
      f"| {'A_90 Xe':>9} {'A_90 Ar':>9}")
print(L)
for fm in ("Fnest", "F1"):
    for t in THR:
        Xr, Ar_ = sx[(t, fm)], sa[(t, fm)]
        print(f"{fm:>6} {t:>6} | {Xr['R_tot']:10.4f} {Ar_['R_tot']:10.4f} "
              f"| {Ar_['R_tot']/Xr['R_tot']:7.1f} | {Xr['A90'][0]:9.4f} {Ar_['A90'][0]:9.4f}")
    print(L)

# Z+2N y Q_W por blanco (de constants.f90; se usan solo en el chequeo analitico)
NSI_GEOM = {"Xe": dict(ZN=208, QW=37.27), "Ar": dict(ZN=62, QW=10.59)}
XLBL, _ = load_nsi_labels("Xe")            # el plano es el mismo para Xe y Ar
DIAG_EE = (r"\varepsilon_{ee}" in XLBL) and (r"^{dV}" in XLBL)

print("\n  CHEQUEO del alcance NSI (por que la proyeccion cubre ~99.8%):")
for tag, S in (("Xe", sx), ("Ar", sa)):
    N = S[(1, "Fnest")]["Ntot192"]
    tol = np.sqrt(4.605 / N)                    # |1 - A_amp| admitido
    print(f"   {tag}:  N_Asimov = {N:6.0f}  ->  tasa CEvNS admitida = SM +/- {100*tol:4.1f} %")
    if DIAG_EE:
        ZN, QW = NSI_GEOM[tag]["ZN"], NSI_GEOM[tag]["QW"]
        r_eps = QW / ZN                         # radio del anillo ciego en eps
        d_eps = QW * tol / (2.0 * ZN)           # semiancho del anillo en eps
        frac_an = 2 * np.pi * r_eps * (2 * d_eps) / 4.0   # area anillo / area caja(=4)
        print(f"        circulo ciego en eps: radio {r_eps:.3f} (grilla: corte eje/2 = "
              f"{blind_info[tag]/2:.3f}), semiancho {d_eps:.4f}  ->  area/caja ~ {frac_an:.4f}")
        print(f"        dentro del alcance (analitico) ~ {100*frac_an:.1f} % ; "
              f"(Fortran, grid 1000x1000) {100-info[tag]:.2f} %")
    else:
        print(f"        (chequeo analitico cerrado solo para el plano ee-diagonal dV; "
              f"aqui vale la grilla)")
    print(f"        => fuera del alcance de la PROYECCION: {info[tag]:.1f} % de la caja |eps|<=1")
print(f"\n  REFERENCIA con datos reales (Xe, rama chi2.f90): A_90 ~ {XE_REAL_A90:.0f} xSM")
print(f"   -> como A_amp max en |eps|<=1 es ~74 < 111, el analisis REAL de Xe NO")
print(f"      distingue del SM NINGUN punto de la caja. El alcance real ~ 0.")

# =====================================================================
# panel resumen 2x2 (opcional, mismo material)
# =====================================================================
fig, AX = plt.subplots(2, 2, figsize=(9.6, 7.6))
(a, b), (c, d) = AX
a.step(ne, rx, where="mid", color="black", lw=1.5)
a.step(ne, ra, where="mid", color="0.45", ls="--", lw=1.5)
a.plot(ne, rx, "o", mfc="white", mec="black", ms=4, ls="none", label="Xe")
a.plot(ne, ra, "s", mfc="0.45", mec="0.45", ms=3.5, ls="none", label="Ar")
a.set_yscale("log"); a.set_xlim(0.4, 21); a.set_ylim(3e-6, 40)
a.set_xlabel(r"$N_e$"); a.set_ylabel(r"$R(N_e)$ [ev/(kg·día)]")
a.set_title("(a) Espectro de ionización"); a.legend()
b.plot(THR, fxv, "o-", color="black", mfc="white", ms=5)
b.plot(THR, fav, "s--", color="0.45", ms=4.5)
b.set_yscale("log"); b.set_xticks(THR); b.set_ylim(0.3, 260)
b.set_xlabel(r"$N_e \geq$"); b.set_ylabel(r"señal retenida [%]")
b.set_title("(b) Coste de subir el umbral"); b.legend(["Xe", "Ar"])
c.plot(expo, sx[(1, "Fnest")]["A90"], "o-", color="black", mfc="white", ms=4.5, label="Xe ideal")
c.plot(expo, sa[(1, "Fnest")]["A90"], "s--", color="0.45", ms=4, label="Ar ideal")
c.axhline(XE_REAL_A90, color="black", lw=1.1, ls=(0, (6, 3)))
c.text(expo[0]*1.1, XE_REAL_A90*0.6, f"Xe real $\\approx$ {XE_REAL_A90:.0f}", fontsize=8)
c.set_xscale("log"); c.set_yscale("log"); c.set_ylim(0.95, 260)
c.set_xlabel(r"exposición [kg·día]"); c.set_ylabel(r"$A_{90}$ [$\times$ SM]")
c.set_title("(c) Sensibilidad proyectada vs exposición"); c.legend()
d.plot(THR, rn, "D-", color="black", mfc="white", ms=5, label="F de NEST")
d.plot(THR, r1, "D--", color="0.5", ms=4.5, label="F=1")
d.axhline(COH_ONLY, color="black", lw=0.9, ls=":")
d.text(2.5, COH_ONLY*1.15, "coherencia sola $\\to$ Xe $\\times$3.8", fontsize=7.5, ha="center")
d.set_yscale("log"); d.set_xticks(THR); d.set_ylim(0.7, 400)
d.set_xlabel(r"$N_e \geq$"); d.set_ylabel(r"eventos Ar / Xe")
d.set_title("(d) Figura de mérito intrínseca"); d.legend()
fig.suptitle("Comparación ideal simétrica Xe vs Ar  ·  sensibilidad PROYECTADA (Asimov)  ·  RED-100",
             fontsize=11.5, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.965], h_pad=2.6, w_pad=2.2)
fig.savefig(f"{BASE}/fig_ideal_XeAr.png")
plt.close(fig)
print(f"\n  {BASE}/fig_ideal_XeAr.png  (resumen 2x2)")

# =====================================================================
# 6) Ar: A_90 vs exposicion en los 3 escenarios de fondo de 39Ar
#    (RED-100 SV con fondo simulado; salida de chi2_bkg_nest.f90)
# =====================================================================
try:
    bk = {}
    for ln in open(f"{BASE}/sensib_bkg_Ar.dat"):
        if ln.lstrip().startswith("#") or not ln.strip():
            continue
        p = ln.split()
        bk.setdefault(p[0], []).append((float(p[1]), float(p[2])))
    fig, ax, cx = new_fig(w=7.8, h=7.6, ratio=(2.3, 1.2))
    sty = {"sin_fondo": ("black", "-", "o", "sin fondo (respuesta intrínseca)"),
           "Ar39_UAr": ("0.35", "--", "s", r"$^{39}$Ar UAr ($7{,}3\times10^{-4}$ Bq/kg)"),
           "Ar39_atmosferico": ("0.6", ":", "^", r"$^{39}$Ar atmosférico (1 Bq/kg)")}
    for key, (col, ls, mk, lab) in sty.items():
        if key not in bk:
            continue
        xs = [e for e, _ in bk[key]]
        ys = [(a - 1.0) for _, a in bk[key]]   # A_90 - 1 (log util)
        ax.plot(xs, ys, ls=ls, marker=mk, color=col, ms=5, mfc="white", label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"exposición reactor ON   [kg$\cdot$día]")
    ax.set_ylabel(r"$A_{90}-1$   [exceso sobre el SM,  90 % C.L.]")
    ax.set_title(r"Ar — sensibilidad esperada (RED-100 §V) con fondo de $^{39}$Ar simulado")
    ax.legend(loc="upper right", fontsize=8.2)
    ax.text(0.03, 0.06,
            r"referencia: Xe con DATOS reales  $A_{90}\!-\!1 \approx 110$"
            "\n" r"(rama chi2.f90, fuera de escala)",
            transform=ax.transAxes, fontsize=8, va="bottom")
    put_caption(cx,
        "Qué muestra  /  Lectura",
        "Método RED-100 §V: señal CEvNS SIMULADA + fondo + Asimov ($\\Delta\\chi^2=2{,}706$).\n"
        "El fondo de $^{39}$Ar se SIMULA (isótopo conocido, $Q_\\beta=565$ keV, actividad\n"
        "publicada) — no se inventa. Escenario 0 = sin fondo (respuesta intrínseca del\n"
        "blanco). El $^{39}$Ar por solape espectral puro apenas mueve $A_{90}$: llegar a\n"
        "$N_e\\!\\leq\\!5$ pide $E_{er}\\!\\lesssim\\!0{,}1$ keV, la cola extrema del espectro $\\beta$.\n"
        "El fondo de apilamiento de electrón único (ref.[46] lo deja SIN RESOLVER)\n"
        "entra como escenario aparte ($s_{\\rm SE}$), no como número: ver README_Ar.txt.\n"
        "NO se reproduce §VI (reactor ON): necesita un detector y un reactor.")
    save(fig, "fig_ideal_6_fondo_Ar.png")
except FileNotFoundError:
    print("  (sensib_bkg_Ar.dat no encontrado; corre chi2_bkg_nest primero)")
