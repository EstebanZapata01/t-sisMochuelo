# RED-100: sensibilidad a NSI y a sin²θ_W, Xe vs. Ar

Código de la tesis de grado de Esteban Zapata (Universidad de Pamplona):
sensibilidad de un detector estilo RED-100 (arXiv:2411.18641, CEνNS de
antineutrinos de reactor) a Interacciones No Estándar y al ángulo de mezcla
débil, comparando un blanco de xenón (validado contra el dato real 2024 de
RED-100) con uno de argón (proyección, ref. [46] + ReD 2025).

**La referencia técnica completa** — toda la física implementada, las
fórmulas, la metodología estadística y la discusión de cada resultado — está
en [`doc/metodologia.tex`](doc/metodologia.tex). Este README solo explica
cómo está organizado el código y cómo correrlo; no repite esa discusión.

## Estructura

```
FORTRAN90/
  N_EventosCEvNS/          rama Ge/CONUS+ original (sin NSI)
  N_EventosCEvNS_NSI/      rama Ge/CONUS+ con NSI 2D y barrido de sin²θ_W
  N_EventosCEvNS_NSIXe/    pipeline RED-100 (xenón): NEST binomial, PE, ideal
  N_EventosCEvNS_NSIAr/    pipeline RED-100 (argón): clon del de Xe + fondo ³⁹Ar
  chi2_nsi_generic/        motor NSI/sin²θ_W genérico (reusado por Xe/Ar/CONUS+)
python/                    post-proceso: figuras, tablas, validaciones, NEST
doc/                       doc/metodologia.tex (fuente) y figuras del pipeline
datos/                     entradas (tablas NEST, datos digitalizados) y
                           salidas .dat/.png/.tex del pipeline (la mayoría se
                           regenera; ver .gitignore)
```

Cada carpeta de `FORTRAN90/` es independiente (mismos módulos compilados por
separado en cada una); no hay un `Makefile` central porque el orden y las
variables de entorno de cada binario están documentados junto a su código —
ver también `FORTRAN90/N_EventosCEvNS_NSIAr/README_Ar.txt` para el registro
detallado de la rama de argón (parámetros, validación contra ref. [46],
resultados de referencia).

## Reproducir todo de un tirón

```bash
./correr_todo.sh          # compila, corre el pipeline y genera las figuras/tablas
./correr_todo.sh figuras  # solo regenera figuras/tablas (asume binarios ya corridos)
```

Al final verifica que estén todas las figuras/tablas citadas en
`doc/metodologia.tex`. Dos cosas a tener en cuenta:

- `python/nest.py` y `python/nest_Ar.py` muestrean NEST por Monte Carlo con
  **semilla fija** (`nestpy.RandomGen.rndm().set_seed(...).lock_seed()`,
  20260910 para Xe y 20260911 para Ar — `np.random.seed()` no sirve porque
  el muestreo es interno a nestpy/C++). Verificado bit-idéntico en corridas
  repetidas. Antes de fijarla, cada corrida completa daba ~0,1–4 % de ruido
  estadístico respecto a la corrida anterior (documentado y discutido con
  la sesión de `Proyecto_LaTeX/` el 2026-09-16); ya no debería ocurrir.
- El script **no** regenera `datos/chi2_sin2theta.dat`: el binario que lo
  produce (`FORTRAN90/N_EventosCEvNS/chi2.f90`, un archivo que este proyecto
  trata como intocable) tiene una variable local (`log_s`) que `gfortran
  -Wuninitialized` marca como posiblemente sin inicializar, y recompilarlo
  desde cero puede dar un $\chi^2(s)$ degenerado (plano) en vez del perfil
  real — ya pasó una vez en esta sesión. El `.dat` se deja versionado, ya
  validado (χ²_min=7,349, s_best=0,248). Si algún día se autoriza tocar ese
  archivo, ahí está la pista para arreglarlo.

## Compilar a mano (gfortran, sin Makefile)

El orden de los módulos en la línea de compilación es obligatorio (cada uno
depende del anterior). Ejemplo para xenón:

```bash
cd FORTRAN90/N_EventosCEvNS_NSIXe
MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal        $MODS chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o chi2red100_nest   $MODS chi2red100_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest   $MODS mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest       $MODS red100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100PE          $MODS mod_detector.f90 red100PE.f90
gfortran -O2 -ffree-line-length-none -o chi2              constants.f90 chi2.f90   # ajuste ON-OFF real (dato 2024)
```

Para argón (mismo patrón, sin `mod_detector.f90`/`red100PE.f90`: el argón no
tiene rama en fotoelectrones, ver `doc/metodologia.tex` §"Por qué el análisis
de argón no pasa por el espectro en fotoelectrones"):

```bash
cd FORTRAN90/N_EventosCEvNS_NSIAr
MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal   $MODS chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o chi2_bkg     $MODS chi2_bkg_nest.f90     # §V con fondo de ³⁹Ar
gfortran -O2 -ffree-line-length-none -o mainred100_nest $MODS mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest  $MODS red100_nest.f90
```

`chi2_nsi_generic/`, `N_EventosCEvNS/` y `N_EventosCEvNS_NSI/` (Ge/CONUS+) se
compilan igual, un binario por programa principal; los módulos compartidos
(`constants.f90`, `flux.f90`, `xsections.f90`, `quenching.f90`,
`resolution.f90`) se listan antes del programa que los usa.

Los binarios, `.mod` y `.o` no se versionan (ver `.gitignore`): son
específicos de cada compilador/máquina y se regeneran con lo de arriba.

## Correr el pipeline y regenerar las figuras a mano (paso a paso)

`./correr_todo.sh` hace justo esto; se detalla por si se quiere correr solo
una parte:

1. Compilar los binarios de la sección anterior.
2. Generar las tablas de NEST que consume Fortran:
   `python python/nest.py` (Xe), `python python/nest_Ar.py` (Ar).
3. Correr los binarios (cada uno escribe sus `.dat` en `datos/`); el orden y
   qué produce cada uno está documentado en el encabezado de cada script
   `.py` de `python/` y en `FORTRAN90/N_EventosCEvNS_NSIAr/README_Ar.txt`.
4. Correr los scripts de `python/` que grafican (todos con
   `MPLBACKEND=Agg python3 python/<script>.py`; cada uno documenta en su
   docstring qué `.dat` lee y qué figura/tabla en `datos/` escribe). El
   estilo visual único de todas las figuras vive en `python/estilo_tesis.py`.

El Apéndice de `doc/metodologia.tex` ("programa → qué produce") es la tabla
de referencia completa entrada→script→salida.

## Convenciones de la física implementada (resumen; detalle en el .tex)

- Todo parámetro numérico se rastrea a una fuente publicada (el paper de
  RED-100, arXiv:2411.18641; ref. [46] + ReD 2025 para argón) o a la salida
  directa del propio código — nada asumido sin citar.
- El fondo de ³⁹Ar (rama de argón) se trata como una constante conocida, no
  como un parámetro libre de la verosimilitud.
- La incertidumbre se explora con corridas de reemplazo discretas
  (`python/sensitivity_scan.py`), al estilo del propio paper de RED-100, no
  con un nuisance continuo.
