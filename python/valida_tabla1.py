#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verificacion del estadistico chi2 con datos ajenos (NO es un resultado propio): reproduce el
limite publicado de RED-100 (SM2018) usando SOLO lo que publica la colaboracion: los residuos
ON-OFF de sus tres histogramas (energia, duracion, radio^2) y su propia prediccion de senal
(barras del limite al 90%, R_i = barra_i/63). Observado: D_i = residuo; Asimov: D_i = R_i.
Usa el mismo binario del ajuste (chi2_nsi_generic).

Entrada: datos/red100_residuo_ONOFF_fig8.csv (python/extrae_pdf_red100.py)
Salidas: datos/tabla_validacion_red100.tex, datos/fig_validacion_red100.png (Delta chi^2(A) reproducido y publicado)
"""
import os
import re
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estilo_tesis import aplicar, C_XE, GRIS, NEGRO, FIG15
aplicar(grande=True)

ROOT = os.path.abspath("..")
ENGINE = f"{ROOT}/FORTRAN90/chi2_nsi_generic"
A90_PUB, ESP_PUB, ABEST_PUB = 63.0, 58.0, 6.0      # Tabla I (SM2018); A* inferido: A90_obs - (A90_esp - 1)

f = [l.strip().split(",") for l in open(f"{ROOT}/datos/red100_residuo_ONOFF_fig8.csv") if not l.startswith(("#", "panel"))]
pan = np.array([r[0] for r in f]); v = np.array([[float(x) for x in r[1:]] for r in f])
D, sig, R = v[:, 1], v[:, 2], v[:, 3] / A90_PUB


def a90_motor(datos, sig, R, tmp, nombre):
    """A_90 (Dchi2 = 2.706) del binario del ajuste con residuos 'datos', sigma y prediccion R."""
    ent = f"{tmp}/{nombre}.dat"
    with open(ent, "w") as g:
        g.write("# Z  N  sigma_alpha  n_bins  ipar\n")
        g.write(f"{54.0:8.2f}{77.3879:10.3f}{-1.0:8.4f}{len(R):6d}{5:6d}\n# bin  dN  sigma  R\n")
        for i, (d, s, r) in enumerate(zip(datos, sig, R), 1):
            g.write(f"{i:5d}{d:16.7e}{s:16.7e}{r:16.7e}\n")
    subprocess.run(["./chi2_nsi_generic", ent, f"{tmp}/{nombre}"], cwd=ENGINE, check=True, capture_output=True)
    return float(re.search(r"A_90\s*=\s*([0-9.]+)", open(f"{tmp}/{nombre}_resumen.txt").read()).group(1))


with tempfile.TemporaryDirectory() as tmp:
    obs, esp = a90_motor(D, sig, R, tmp, "obs"), a90_motor(R, sig, R, tmp, "asimov")
    E = pan == "E"                                   # solo el histograma de energia, con la senal publicada
    obs_E, esp_E = a90_motor(D[E], sig[E], R[E], tmp, "obsE"), a90_motor(R[E], sig[E], R[E], tmp, "asimovE")
S1, S2 = (D * R / sig ** 2).sum(), (R ** 2 / sig ** 2).sum()
S2E = (R[pan == "E"] ** 2 / sig[pan == "E"] ** 2).sum()

print("  Tabla I (SM2018)      publicado   reproducido")
print(f"  A_90 observado        {A90_PUB:8.0f}   {obs:9.1f}")
print(f"  A_90 esperado         {ESP_PUB:8.0f}   {esp:9.1f}")
print(f"  A* = S1/S2            {ABEST_PUB:8.0f}   {S1 / S2:9.1f}   (publicado: inferido de la Tabla I)")
print(f"  S2(3 hist)/S2(energia)          {S2 / S2E:9.2f}")
print(f"  solo energia, senal publicada:  A_90 obs {obs_E:.1f}  esp {esp_E:.1f}  (factor 3 hist. en A_90-1: {(esp_E - 1) / (esp - 1):.2f})")

c = lambda x, d=1: f"{x:.{d}f}".replace(".", "{,}")          # coma decimal
with open(f"{ROOT}/datos/tabla_validacion_red100.tex", "w", encoding="utf-8") as g:
    g.write("% generado por python/valida_tabla1.py\n\\begin{tabular}{@{}lcc@{}}\n\\toprule\n"
            "& publicado & reproducido \\\\\n\\midrule\n")
    g.write(f"$A_{{90}}$ observado [$\\times$SM] & ${c(A90_PUB, 0)}$ & ${c(obs)}$ \\\\\n")
    g.write(f"$A_{{90}}$ esperado [$\\times$SM] & ${c(ESP_PUB, 0)}$ & ${c(esp)}$ \\\\\n")
    g.write(f"$A^{{*}}=S_1/S_2$ & $\\approx{c(ABEST_PUB, 0)}$ & ${c(S1 / S2)}$ \\\\\n\\bottomrule\n\\end{{tabular}}\n")
print("  -> datos/tabla_validacion_red100.tex")

# ---- figura: Delta chi^2(A) reproducido frente al publicado (curvas extraidas del PDF)
pub = {}
for l in open(f"{ROOT}/datos/red100_perfil_chi2_fig9.csv"):
    if not l.startswith(("#", "curva")):
        c, a_, d_ = l.strip().split(",")
        pub.setdefault(c, []).append((float(a_), float(d_)))
S3 = (D ** 2 / sig ** 2).sum()
A = np.linspace(0, 100, 1001)
chi2_obs = S3 - 2 * A * S1 + A ** 2 * S2
d_obs = chi2_obs - chi2_obs[np.argmin(np.where(A >= 0, chi2_obs, np.inf))]        # minimo fisico (A >= 0)
d_esp = S2 * (A - 1) ** 2
fig, ax = plt.subplots(figsize=FIG15)
ax.axhline(2.706, color=NEGRO, lw=1.1, label=r"corte al $90\%$ C.L. ($\Delta\chi^{2}=2{,}71$)")
ax.plot(A, d_esp, color=GRIS, ls="-.", label=f"esperado, reproducido: $A_{{90}}={esp:.1f}$".replace(".", "{,}"))
ax.plot(A, d_obs, color=C_XE, label=f"observado, reproducido: $A_{{90}}={obs:.1f}$".replace(".", "{,}"))
for k, ls in (("esperada", ":"), ("observada", ":")):
    q = np.array(pub[k])
    ax.plot(q[:, 0], q[:, 1], color=NEGRO, ls=ls, lw=1.6, label="publicado" if k == "esperada" else None)
for x, col in ((esp, GRIS), (obs, C_XE)):
    ax.vlines(x, 0, 2.706, color=col, ls=":", lw=1.4)
    ax.plot([x], [2.706], "o", mfc="white", mec=col, mew=1.6, ms=8, zorder=5)
ax.set_xlim(0, 100); ax.set_ylim(0, 5)
ax.set_xlabel(r"amplitud $A$  [$\times$SM]")
ax.set_ylabel(r"$\Delta\chi^{2}$")
ax.legend(loc="upper left")
fig.tight_layout(); fig.savefig(f"{ROOT}/datos/fig_validacion_red100.png"); plt.close(fig)
print("  -> datos/fig_validacion_red100.png")
