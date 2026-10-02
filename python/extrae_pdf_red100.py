#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extrae de las figuras vectoriales de RED-100 (arXiv:2411.18641v1) los datos publicados, leyendo la
geometria del PDF (pdftocairo -svg) y calibrando cada eje con sus marcas:
  pag. 7  residuos ON-OFF (energia, duracion, radio^2), sus barras de error y las barras del limite al 90%
  pag. 8  curvas Delta chi^2(A) observada (ON-OFF) y esperada (Asimov)
  pag. 6  senal CEvNS y fondo medido por bin de N_e reconstruido, antes y despues de los cortes

Uso (herramienta puntual; los CSV resultantes van al repo):  python3 extrae_pdf_red100.py [pdf]
Salidas: datos/red100_residuo_ONOFF_fig8.csv  (panel, x, dN, sigma, limite90; panel E = energia [PE], D = duracion [ns], R = radio^2 [mm^2])
         datos/red100_perfil_chi2_fig9.csv    (curva, A, dchi2)
         datos/red100_senal_fondo_fig6.csv    (Ne, senal_antes, senal_despues, fondo_antes, fondo_despues) en cuentas/(kg dia)
"""
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
import numpy as np

ROOT = os.path.abspath("..")
PDF = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/RED100_2411.18641v1.pdf"
PAG = 7
A90_PUB = 63.0                                    # limite publicado (SM2018), con el que se dibujan las barras
PANEL = {"E": (243, 305), "D": (338, 400), "R": (433, 495)}          # rango en y de cada panel (coord. de pagina)
XTICK = {"E": [110 + 10 * i for i in range(9)], "D": [2000, 2500, 3000, 3500, 4000], "R": [0, 5000, 10000, 15000]}
NS = "{http://www.w3.org/2000/svg}"


def nums(s):
    return [float(v) for v in re.findall(r"-?\d+\.?\d*(?:e-?\d+)?", s)]


def lee_paths(svg):
    """(estilo, puntos en coordenadas de pagina) de cada <path> fuera de <defs>."""
    out = []

    def rec(e):
        for c in e:
            t = c.tag.replace(NS, "")
            if t == "defs":
                continue
            if t == "path":
                v = nums(c.get("d") or "")
                p = np.array(list(zip(v[0::2], v[1::2])))
                tr = c.get("transform")
                if tr and len(p):
                    a, b, cc, d, e0, f = nums(tr)
                    p = np.column_stack([a * p[:, 0] + cc * p[:, 1] + e0, b * p[:, 0] + d * p[:, 1] + f])
                if len(p):
                    out.append((c.get("style") or "", p))
            rec(c)

    rec(ET.parse(svg).getroot())
    return out


def panel_de(y):
    return next((k for k, (lo, hi) in PANEL.items() if lo <= y <= hi), None)


def pagina(n):
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdftocairo", "-svg", "-f", str(n), "-l", str(n), PDF, f"{tmp}/p.svg"], check=True)
        return lee_paths(f"{tmp}/p.svg")


paths = pagina(PAG)

barras = {k: [] for k in PANEL}; err = {k: [] for k in PANEL}; marc = {k: [] for k in PANEL}
xt = {k: [] for k in PANEL}; yt = {k: [] for k in PANEL}
for st, p in paths:
    x0, x1, y0, y1 = p[:, 0].min(), p[:, 0].max(), p[:, 1].min(), p[:, 1].max()
    k = panel_de((y0 + y1) / 2)
    if "49.803162%" in st and k and 3 < x1 - x0 < 60 and 340 < x0 and x1 < 565:      # barras naranja
        barras[k].append((x0, x1, y0, y1))
    elif "stroke-width:1.5" in st and "stroke:rgb(0%,0%,0%)" in st and k and x1 - x0 < 0.05:   # barra de error
        err[k].append((x0, y0, y1))
    elif "fill:rgb(0%,0%,0%)" in st and "fill:none" not in st and k and 2.5 < x1 - x0 < 7 and y1 - y0 < 7 and 340 < x0 and x1 < 565:
        marc[k].append(((x0 + x1) / 2, (y0 + y1) / 2))                                # marcador
    elif "stroke-width:0.8" in st and "stroke:rgb(0%,0%,0%)" in st and len(p) == 2 and 335 < p[:, 0].min() < 565:
        dx, dy = abs(p[0, 0] - p[1, 0]), abs(p[0, 1] - p[1, 1])
        if dy < 0.01 and 1 < dx < 3 and panel_de(p[0, 1]):                            # marca en y
            yt[panel_de(p[0, 1])].append(p[0, 1])
        elif dx < 0.01 and 1 < dy < 3:                                                # marca en x
            for kk, (lo, hi) in PANEL.items():
                if hi < p[:, 1].min() < hi + 8:
                    xt[kk].append(p[0, 0])

filas = []
for k in "EDR":
    tx, ty = np.array(sorted(xt[k])), np.array(sorted(yt[k]))          # ty: marcas 0.2, 0.0, -0.2 (de arriba a abajo)
    ax, bx = np.polyfit(tx, XTICK[k], 1)
    escala, y_cero = (ty[1] - ty[0]) / 0.2, ty[1]                       # pt por (cuentas/kg/dia); y = 0
    mk = sorted(marc[k])
    dx = np.median(np.diff([m[0] for m in mk]))
    mk = [m for m in mk if min(((m[0] - mk[0][0]) / dx) % 1, 1 - ((m[0] - mk[0][0]) / dx) % 1) < 0.15]   # sin el marcador de la leyenda
    bar = [b for b in barras[k] if abs(b[3] - y_cero) < 0.6]                                              # sin el parche de la leyenda
    assert len(mk) == len(bar) == 15, (k, len(mk), len(bar))
    for (mx, my), b in zip(mk, sorted(bar)):
        e = next(v for v in sorted(err[k]) if abs(v[0] - mx) < 0.8)
        filas.append((k, ax * mx + bx, (y_cero - my) / escala, (e[2] - e[1]) / 2 / escala, (y_cero - b[2]) / escala))

cab = (f"# Residuo ON-OFF de RED-100 y limite al 90% C.L. por bin (arXiv:2411.18641v1, pagina {PAG}, figura de residuos).\n"
       "# Metodo: pdftocairo -svg sobre la pagina; geometria vectorial exacta (marcadores, barras de error, barras naranja);\n"
       "# ejes calibrados con sus marcas (E: 110-190 PE; D: 2000-4000 ns; R: 0-15000 mm^2; y: 0.2, 0.0, -0.2 cuentas/(kg dia)).\n"
       f"# limite90 = A_90 * R_i con A_90 = {A90_PUB:g} (SM2018, Tabla I del paper): R_i = limite90/{A90_PUB:g}.\n"
       "# Paneles: E energia corregida [PE], D duracion [ns], R radio^2 [mm^2]. sigma = semilongitud de la barra de error.\n"
       "# Generado por python/extrae_pdf_red100.py\n")
with open(f"{ROOT}/datos/red100_residuo_ONOFF_fig8.csv", "w", encoding="utf-8") as f:
    f.write(cab + "panel,x,dN,sigma,limite90\n")
    for r in filas:
        f.write(f"{r[0]},{r[1]:.4f},{r[2]:.6f},{r[3]:.6f},{r[4]:.6f}\n")
print(f"  -> datos/red100_residuo_ONOFF_fig8.csv ({len(filas)} filas)")


# ------------------------------------------------------------------ curvas Delta chi^2(A) (pag. 8)
def marcas(paths, eje, rango, fija):
    """Posiciones de las marcas (trazo negro de 0.8, largo 1-3) de un eje dentro de 'rango' con la otra coordenada en 'fija'."""
    out = []
    for st, p in paths:
        if "stroke-width:0.8" in st and "stroke:rgb(0%,0%,0%)" in st and len(p) == 2:
            q = p[0, eje]
            if rango[0] < q < rango[1] and fija[0] < p[:, 1 - eje].min() < fija[1] and abs(p[0, 1 - eje] - p[1, 1 - eje]) > 1 and abs(p[0, eje] - p[1, eje]) < 0.01:
                out.append(round(q, 3))
    return sorted(set(out))


p8 = pagina(8)
tx = marcas(p8, 0, (75, 265), (140, 152))[:5]                          # A = 0, 20, 40, 60, 80
ty = [q for q in sorted(set(round(p[0, 1], 3) for st, p in p8 if "stroke-width:0.8" in st and "stroke:rgb(0%,0%,0%)" in st and len(p) == 2
                            and abs(p[0, 1] - p[1, 1]) < 0.01 and 70 < p[0, 0] < 82 and 100 < p[0, 1] < 150))]   # marcas 0, 2 (y en pagina)
y0, y2 = 144.621, 109.137                                              # Dchi2 = 0 y 2 (marcas de los ejes de la figura)
assert abs(ty[-1] - y0) < 0.01 and any(abs(q - y2) < 0.01 for q in ty)
with open(f"{ROOT}/datos/red100_perfil_chi2_fig9.csv", "w", encoding="utf-8") as f:
    f.write("# Curvas Delta chi^2(A) de RED-100 (arXiv:2411.18641v1, pagina 8): observada (ON-OFF) y esperada (Asimov), SM2018.\n"
            "# Metodo: pdftocairo -svg; vertices de las dos curvas; ejes calibrados con sus marcas (A: 0-80 x SM; Dchi2: 0 y 2).\n"
            "# Generado por python/extrae_pdf_red100.py\ncurva,A,dchi2\n")
    for nombre, clave in (("observada", "58.03833%"), ("esperada", "100%,49.803162%")):
        p = next(pp for st, pp in p8 if len(pp) > 20 and clave in st and "stroke-width" in st)
        for x, y in sorted(zip(80.0 * (p[:, 0] - tx[0]) / (tx[-1] - tx[0]), 2.0 * (y0 - p[:, 1]) / (y0 - y2))):
            f.write(f"{nombre},{x:.4f},{y:.5f}\n")

# ------------------------------------------------------------------ senal y fondo por N_e (pag. 6)
p6 = pagina(6)


def cuadro(paths, clave, ylo, yhi):
    """Centros (x, y) de los marcadores de un color dentro de la banda vertical [ylo, yhi]."""
    out = []
    for st, p in paths:
        if clave in st and "fill:none" not in st and 2.5 < np.ptp(p[:, 0]) < 6 and 365 < p[:, 0].mean() < 560:
            cy = (p[:, 1].min() + p[:, 1].max()) / 2
            if ylo < cy < yhi:
                out.append((p[:, 0].mean(), cy))
    return out


def decadas(ylo, yhi):
    return sorted(set(round(p[0, 1], 3) for st, p in p6 if "stroke:rgb(0%,0%,0%)" in st and len(p) == 2 and abs(p[0, 1] - p[1, 1]) < 0.01
                      and 360 < p[:, 0].min() < 372 and ylo < p[0, 1] < yhi and abs(abs(p[0, 0] - p[1, 0]) - 1.75) < 0.2))


xcol = sorted(set(round(p[0, 0], 1) for st, p in p6 if "stroke-width:0.8" in st and "stroke:rgb(0%,0%,0%)" in st and len(p) == 2
                  and abs(p[0, 0] - p[1, 0]) < 0.01 and 262 < p[:, 1].min() < 270 and 365 < p[0, 0] < 560))[:4]      # N_e = 4..7


def serie(paths, clave, dec, ylo, yhi, v0):
    y0, d = dec[0], (dec[-1] - dec[0]) / (len(dec) - 1)
    m = cuadro(paths, clave, ylo, yhi)
    return [v0 * 10 ** (-(next(c for x, c in m if abs(x - xc) < 0.6) - y0) / d) for xc in xcol]


dsen, dfon = decadas(180, 265), decadas(100, 182)                       # marcas mayores: 1e-1..1e-4 y 1e3..1e0
sa = serie(p6, "fill:rgb(0%,0%,0%)", dsen, 185, 262, 1e-1)               # senal antes de cortes (circulos negros)
sd = serie(p6, "49.803162%", dsen, 185, 262, 1e-1)                       # senal despues (triangulos naranja)
fa = serie(p6, "fill:rgb(0%,0%,0%)", dfon, 100, 182, 1e3)                # fondo antes (cuadrados negros)
fd = serie(p6, "83.920288%", dfon, 100, 182, 1e3)                        # fondo despues (triangulos rojos)
with open(f"{ROOT}/datos/red100_senal_fondo_fig6.csv", "w", encoding="utf-8") as f:
    f.write("# Senal CEvNS simulada (SM2018) y fondo medido (reactor apagado) por bin de N_e reconstruido, antes y despues de los cortes,\n"
            "# en cuentas/(kg dia) (arXiv:2411.18641v1, pagina 6). Metodo: pdftocairo -svg; ejes calibrados con las marcas mayores de los ejes.\n"
            "# Generado por python/extrae_pdf_red100.py\nNe,senal_antes,senal_despues,fondo_antes,fondo_despues\n")
    for i, n in enumerate(range(4, 8)):
        f.write(f"{n},{sa[i]:.6e},{sd[i]:.6e},{fa[i]:.6e},{fd[i]:.6e}\n")
print("  -> datos/red100_perfil_chi2_fig9.csv, datos/red100_senal_fondo_fig6.csv")
