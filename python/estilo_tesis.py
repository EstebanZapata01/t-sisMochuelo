#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
estilo_tesis.py -- estilo unico para todas las figuras de la metodologia.

Paleta OPACA (baja saturacion, apta para impresion y legible en escala de
grises: los tonos difieren tambien en luminosidad). Tipografia serif, ejes
finos, marcas hacia adentro, sin rejilla pesada ni marcos de leyenda.

Uso en cada script de figura, justo despues de importar pyplot:

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from estilo_tesis import aplicar, C_XE, C_AR   # etc.
    aplicar()
"""
import matplotlib as _mpl
from cycler import cycler as _cycler

# --------------------------------------------------------------- paleta opaca
NEGRO   = "#1a1a1a"
GRIS    = "#8a8a8a"
GRIS_CL = "#c7c7c7"

AZUL    = "#33546e"   # azul acero apagado   (Xe / serie 1)   ~ reemplaza #1f77b4
NARANJA = "#a86a43"   # terracota apagado    (Ar / serie 2)   ~ reemplaza #E87722
VERDE   = "#5c7053"   # verde salvia apagado (serie 3)        ~ reemplaza #2ca02c
ROJO    = "#8f4444"   # ladrillo apagado     (serie 4)        ~ reemplaza #d62728
MORADO  = "#6a5a78"   # ciruela apagado      (serie 5)        ~ reemplaza #9467BD
MARRON  = "#7a5c4a"   # marron apagado       (serie 6)
NAVY    = "#2b3f52"   # azul noche           (acento oscuro)
ARENA   = "#c89a6a"   # arena apagado        (serie clara aux)
TEAL    = "#4f7068"   # verde azulado apagado

# alias por nombre de dominio (muchos scripts los usan asi)
C_XE, C_AR, C_GE = AZUL, NARANJA, VERDE
C_SM = NEGRO

CICLO = [AZUL, NARANJA, VERDE, ROJO, MORADO, MARRON]

# mapa secuencial para rellenos de contornos / mapas 2D
CMAP  = "Greys"


def aplicar():
    """Fija los rcParams del estilo de la tesis (idempotente)."""
    _mpl.rcParams.update({
        # tipografia
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "font.size": 10,
        "axes.titlesize": 10.5,
        "axes.labelsize": 10,
        "legend.fontsize": 8.6,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        # ejes y marcas
        "axes.linewidth": 0.8,
        "axes.edgecolor": NEGRO,
        "axes.labelcolor": NEGRO,
        "text.color": NEGRO,
        "xtick.color": NEGRO,
        "ytick.color": NEGRO,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "lines.linewidth": 1.6,
        "lines.markersize": 5,
        # rejilla discreta (si el script la activa)
        "grid.color": GRIS_CL,
        "grid.linewidth": 0.4,
        "grid.linestyle": ":",
        "axes.grid": False,
        # leyenda sobria
        "legend.frameon": False,
        # color por defecto de las series
        "axes.prop_cycle": _cycler(color=CICLO),
        # figura
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "image.cmap": CMAP,
    })
