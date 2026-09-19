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
| `sim_extraidos_x_effROI` | `sim_extraidos` × `effROI` (lo que entra al ajuste), solo Ne = 4–7 |
| `red100_senal_antes_cortes`, `red100_senal_despues_cortes` | señal CEνNS simulada por el paper antes/después de sus cortes (Fig. 6 abajo), digitalizada, Ne = 4–7 |
| `effROI` | retención por bin usada en el ajuste (0,138; 0,330; 0,599; 0,719) |

Los valores de RED-100 los digitalizó el autor con WebPlotDigitizer sobre las
figuras (ejes logarítmicos). Dos pasadas independientes sobre la Fig. 6 difieren
como máximo 3 % bin a bin. `effROI` viene de una pasada anterior: el cociente
después/antes de las columnas actuales da 0,137; 0,327; 0,611; 0,737.

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

`sim_creados`/`sim_extraidos` no llevan eficiencia de selección ni tiempo vivo.
Al triplicar los nodos en T (6000), los valores cambian ≲ 0,1 %.

## Antes de comparar

- **La Fig. 3 y la Fig. 6 (antes de cortes) del paper no son la misma curva**:
  en N_e = 4 difieren casi ×2 y las razones entre bins consecutivos son otras
  (~7 contra ~1,2 entre N_e = 4 y 5). Según el texto del paper, la de la
  Fig. 6 ya pasó por su reconstrucción de posición, duración y energía
  corregida; la de la Fig. 3 no. Esa reconstrucción no se puede reproducir con
  lo publicado.
- La simulación sigue la forma de la Fig. 3 con un factor ~1,8–2,1 en Ne = 4–7.
  Se atribuye (sin verificar) a que el flujo híbrido es más duro que el SM2018
  del paper.
- No uses Ne = 0 para validar: la simulación solo cuenta retrocesos con
  T ≥ 0,2 keV y no se sabe qué incluye el paper en ese bin.
- Después de los cortes la simulación queda ×3,8 arriba en Ne = 4 y ×0,6–0,7 en
  Ne = 5–7 respecto a la señal del paper.

## No incluido

Los datos de fondo (paneles superiores de la Fig. 6): existe una digitalización
preliminar, pero sus puntos de señal no coinciden con la pasada posterior, así
que no se garantiza y no se incluye.
