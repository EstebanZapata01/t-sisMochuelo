#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
estilo_tesis.py -- estilo unico para TODAS las figuras de la metodologia.

Principios (aplicarlos, no solo importar la paleta):
  * poco ruido: sin rejilla pesada, sin marcos de leyenda, spines finos,
    marcas hacia adentro y cortas.
  * pocos marcadores: lineas limpias; marcador solo donde hay puntos
    discretos reales (datos, no la curva).
  * poco texto DENTRO de la grafica: el titulo es opcional y corto, la
    interpretacion va en el pie de figura del .tex. Anotaciones = flecha
    fina + 2-4 palabras (helper `nota`), nunca un recuadro con parrafo.
  * el "zoom" a un detalle se hace con un inset y lineas conectoras
    (helper `zoom`), no con un segundo panel suelto.
  * paleta OPACA, apta para impresion y escala de grises (los tonos
    difieren tambien en luminosidad).

Uso:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from estilo_tesis import aplicar, C_XE, C_AR, nota, zoom, sombra
    aplicar()
"""
import matplotlib as _mpl
from cycler import cycler as _cycler

# --------------------------------------------------------------- paleta opaca
NEGRO   = "#1a1a1a"
GRIS    = "#8a8a8a"
GRIS_CL = "#c9c9c9"

AZUL    = "#33546e"   # azul acero apagado   -> Xe / serie 1
NARANJA = "#a86a43"   # terracota apagado    -> Ar / serie 2
VERDE   = "#5c7053"   # verde salvia apagado -> Ge / serie 3
ROJO    = "#8f4444"   # ladrillo apagado     -> serie 4
MORADO  = "#6a5a78"   # ciruela apagado      -> serie 5
MARRON  = "#7a5c4a"   # marron apagado       -> serie 6
NAVY    = "#2b3f52"
ARENA   = "#c89a6a"
TEAL    = "#4f7068"

C_XE, C_AR, C_GE = AZUL, NARANJA, VERDE
C_SM = NEGRO
CICLO = [AZUL, NARANJA, VERDE, ROJO, MORADO, MARRON]
CMAP  = "Greys"


def aplicar():
    """Fija los rcParams de la tesis (idempotente)."""
    _mpl.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "font.size": 10.0,
        "axes.titlesize": 10.5,
        "axes.titlepad": 8.0,
        "axes.labelsize": 10.0,
        "legend.fontsize": 8.4,
        "xtick.labelsize": 8.8,
        "ytick.labelsize": 8.8,
        # ejes: finos, sobrios, contenido por encima de la rejilla
        "axes.linewidth": 0.8,
        "axes.edgecolor": NEGRO,
        "axes.labelcolor": NEGRO,
        "axes.axisbelow": True,
        "text.color": NEGRO,
        "xtick.color": NEGRO, "ytick.color": NEGRO,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        # trazos
        "lines.linewidth": 1.7,
        "lines.markersize": 4.0,
        "lines.markeredgewidth": 0.0,
        "patch.linewidth": 0.0,
        # rejilla (solo si el script la pide con ax.grid())
        "grid.color": GRIS_CL, "grid.linewidth": 0.4, "grid.linestyle": ":",
        "axes.grid": False,
        # leyenda
        "legend.frameon": False,
        "legend.handlelength": 1.7,
        "legend.handletextpad": 0.6,
        "legend.borderaxespad": 0.5,
        "legend.labelspacing": 0.35,
        # series
        "axes.prop_cycle": _cycler(color=CICLO),
        # figura / export
        "figure.facecolor": "white",
        "figure.dpi": 120,
        "savefig.facecolor": "white",
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
        "image.cmap": CMAP,
    })


def limpiar(ax):
    """Quita los spines superior y derecho (look mas ligero, opcional)."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(top=False, right=False)


def nota(ax, xy, text, xytext, color=NEGRO, fs=8.2, **kw):
    """Anotacion estandar: flecha fina + 2-4 palabras, sin recuadro."""
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=fs, color=color,
                arrowprops=dict(arrowstyle="-", lw=0.7, color=color,
                                shrinkA=0, shrinkB=2), **kw)


def sombra(ax, x0, x1, color=GRIS, alpha=0.10, **kw):
    """Banda vertical sutil (ROI, valle de estabilidad, ...)."""
    return ax.axvspan(x0, x1, color=color, alpha=alpha, lw=0, zorder=0, **kw)


def zoom(ax_parent, bounds, xlim, ylim, edge=GRIS, lw=0.8):
    """
    Inset de zoom con lineas conectoras al rectangulo indicado.
    bounds = [x0, y0, w, h] en fraccion de ejes del panel padre.
    Devuelve el eje del inset (ya con el estilo aplicado).
    """
    axin = ax_parent.inset_axes(bounds)
    axin.set_xlim(*xlim)
    axin.set_ylim(*ylim)
    axin.tick_params(labelsize=7, length=2.2, width=0.6)
    for s in axin.spines.values():
        s.set_linewidth(0.7)
        s.set_edgecolor(edge)
    axin.set_facecolor("white")
    ax_parent.indicate_inset_zoom(axin, edgecolor=edge, linewidth=lw, alpha=0.9)
    return axin


def cl_lines(ax, niveles=((1.0, r"$1\sigma$"), (2.706, r"$90\%$"), (3.84, r"$2\sigma$")),
             x_label=1.004):
    """Lineas horizontales de nivel de confianza para paneles Delta chi^2."""
    ymax = ax.get_ylim()[1]
    for y, lab in niveles:
        if y > ymax:
            continue
        ax.axhline(y, color=GRIS, lw=0.7, ls="-" if abs(y - 2.706) < 1e-6 else ":")
        ax.text(x_label, y, lab, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=7.4, color=GRIS)
