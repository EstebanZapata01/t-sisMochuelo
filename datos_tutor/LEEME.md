# Datos de xenón (RED-100) para verificación independiente

Todo sale del pipeline de este repositorio; nada está escrito a mano.
Se regenera con `python/exporta_datos_tutor.py` (o `./correr_todo.sh`).

## Tabla principal: `espectros_Xe.csv`

`Ne` = número de electrones de ionización (0 a 15). Unidades: eventos/(kg·día).
Celda vacía = no disponible.

| columna | qué es |
|---|---|
| `sim_creados`, `sim_extraidos` | simulación propia, antes y después de la extracción (`red100_nest.f90`) |
| `red100_creados`, `red100_extraidos` | espectro de N_e del paper de RED-100 (arXiv:2411.18641, Fig. 3 arriba), digitalizado |
| `sim_reconstruido` | simulación en N_e **reconstruido**: espectro en PE del código (`ionization_spectra_detallado.dat`), PE/27, ventanas de ±0,5 e⁻ recortadas a la ROI 110–189 PE; Ne = 4–7 |
| `sim_ajuste_effROI` | `sim_reconstruido` × `effROI`: predicción de la señal después de cortes; Ne = 4–7 |
| `red100_senal_antes_cortes`, `red100_senal_despues_cortes` | señal CEνNS simulada por el paper antes/después de sus cortes (Fig. 6 abajo), digitalizada, Ne = 4–7 |
| `effROI` | retención por bin de N_e **reconstruido** usada en el ajuste (0,1369; 0,3271; 0,6105; 0,7373) |

Los valores de RED-100 los digitalizó el autor con WebPlotDigitizer sobre las
figuras (ejes logarítmicos). Dos pasadas independientes sobre la Fig. 6 difieren
como máximo 3 % bin a bin. `effROI` es el cociente después/antes de esas mismas columnas.

## Cómo se calcula la simulación (para reproducirla)

1. **Flujo** (`insumos/flujo_nu_hibrido.csv`): Kopeikin para E_ν < 2 MeV
   (tabla, interpolación lineal) y Huber–Mueller para E_ν ≥ 2 MeV (suma sobre
   U235/U238/Pu239/Pu241, fracciones de fisión 0,717/0,068/0,184/0,031),
   normalizado a Φ_total = 1,4×10¹³ ν̄ cm⁻² s⁻¹.
2. **Sección eficaz**: dσ/dT = (G_F² M/π) Q_W² (1 − MT/2E_ν² − T/E_ν), con
   F² = 1. Z = 54, A = 131,293 (peso atómico), N = A − Z = 77,293,
   M = A·931,49410242 MeV, sin²θ_W = 0,23857,
   Q_W = −N/2 + (1 − 4 sin²θ_W) Z/2 = −37,4121, G_F = 1,1663787×10⁻¹¹ MeV⁻²,
   (ħc)² = 3,89379×10⁻²² cm² MeV².
3. **Tasa de retroceso** (`insumos/retroceso_Xe.csv`, columna `dRdT_hibrido`,
   eventos/(kg·día·keV)): ∫ Φ(E_ν) dσ/dT dE_ν sobre E_ν de √(MT/2) a 10 MeV,
   por (N_A/A·1000) átomos/kg y 86400 s/día.
4. **Ionización** (`insumos/nest_Xe_Qy_F.csv`, NEST 2.4.5, 218 V/cm,
   ρ = 2,96 g/cm³, interpolación lineal): ⟨N_e⟩ = T·Q_y(T); p_F = 1 − F(T) con F
   acotado a [0,02; 0,98]; n_F = nint(⟨N_e⟩/p_F) (mínimo 1).
   Creados ~ Bin(n_F, p_F); extraídos ~ Bin(n_F, p_F·EEE), EEE = 0,328.
5. **Espectro en N_e**: suma sobre T ∈ [0,2; 2] keV (trapecio, 2000 nodos) de
   dR/dT · ΔT · P(N_e | T). N_e se trunca en 15. El piso de 0,2 keV es el de la
   tabla NEST; el retroceso máximo físico es 1,64 keV.

6. **Respuesta en PE y N_e reconstruido** (`red100PE.f90`, `mod_detector.f90`): k electrones
   verdaderos → PE ~ Normal(27·k, √k·7,6). N_e reconstruido = PE/27 al entero más cercano; el
   bin j (4..7) es la ventana de PE [(j−½)·27, (j+½)·27] recortada a 110–189 PE. Probabilidad de
   migración P(k→j) = Φ((b_j−27k)/σ_k) − Φ((a_j−27k)/σ_k). `sim_reconstruido`(j) =
   Σ_{k=1..15} R_k · P(k→j) con R_k = `sim_extraidos`; `sim_ajuste_effROI`(j) = `effROI`(j) ·
   `sim_reconstruido`(j). Estas dos columnas salen de `red100PE.f90` (malla de 500 nodos en T
   hasta 3 keV); difieren ≲ 1 % de `sim_extraidos`, que usa 2000 nodos hasta 2 keV.

`sim_creados`/`sim_extraidos` no llevan eficiencia de selección ni tiempo vivo.
Al triplicar los nodos en T (6000), los valores cambian ≲ 0,1 %.

## Antes de comparar

- **Hay dos espacios y no se mezclan.**
  1. *N_e verdadero*: `sim_creados`, `sim_extraidos` frente a `red100_creados`, `red100_extraidos`.
  2. *N_e reconstruido* (PE corregido / 27, la definición del paper): `sim_reconstruido`,
     `sim_ajuste_effROI` frente a `red100_senal_*`. La resolución en PE mueve eventos entre
     bins, así que este espectro es mucho más plano que el verdadero (razón Ne = 4→5 de ~7
     a ~1,2). Comparar un espacio con el otro da conclusiones falsas.
- El agrupamiento en ventanas de ±0,5 e⁻ (cada evento al entero más cercano a PE/27), recortadas
  a 110–189 PE, es un **supuesto**: el paper solo da que 4–7 e⁻ equivalen a 110–189 PE y no
  detalla cómo agrupa. Su respaldo es que la forma resultante coincide con la señal publicada.
- **Espacio reconstruido**: la forma coincide (razones respecto a Ne = 4: 0,80 / 0,137 / 0,016
  la simulación frente a 0,83 / 0,143 / 0,017 RED-100). Cocientes sim / RED-100: 1,33; 1,28;
  1,28; 1,21, iguales antes y después de cortes (`effROI` es el mismo cociente). El factor
  1,2–1,3 no está explicado.
- **Espacio verdadero**: la simulación sigue la forma de los extraídos con un factor ~1,8–2,1
  en Ne = 4–7, **sin explicar**. Los *creados* sí coinciden con los del paper (suma
  Ne = 1–10: 24,2 frente a 24,6), así que el desfase aparece en la extracción. Observación
  sobre los datos digitalizados: los extraídos del paper equivalen a adelgazar sus propios
  creados con p ≈ 0,29 en vez de EEE = 0,328 (≈ 0,89·EEE); no se identificó el origen. No se
  adopta: en el espacio reconstruido ese factor empeora la forma (razones 0,755 / 0,118 /
  0,013 frente a 0,83 / 0,143 / 0,017), y la EEE del paper se midió con calibración gamma.
- No uses Ne = 0. Además, los extraídos digitalizados no conservan el número de eventos
  (Σ extraídos = 22,7 frente a Σ creados = 33,4 con el punto Ne = 0 leído), lo que indica un
  posible error de digitalización en ese punto.

## No incluido

Los datos de fondo (paneles superiores de la Fig. 6): existe una digitalización
preliminar, pero sus puntos de señal no coinciden con la pasada posterior, así
que no se garantiza y no se incluye.
