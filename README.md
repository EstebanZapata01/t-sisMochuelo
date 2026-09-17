# RED-100: sensibilidad a NSI y a sin²θ_W, Xe vs. Ar

Pipeline de la tesis de grado de Esteban Zapata (Universidad de Pamplona).
Compara la sensibilidad de un detector estilo RED-100 (arXiv:2411.18641,
CEνNS de antineutrinos de reactor) a NSI y al ángulo de mezcla débil, para
un blanco de xenón (validado contra el dato real 2024) y uno de argón
(proyección, ref. [46] + ReD 2025).

La física y la metodología completas están en
[`doc/metodologia.tex`](doc/metodologia.tex). Este README solo explica
cómo está organizado el código y cómo correrlo.

## Estructura

```
FORTRAN90/
  N_EventosCEvNS/          Ge/CONUS+ original (sin NSI)
  N_EventosCEvNS_NSI/      Ge/CONUS+ con NSI 2D y barrido de sin²θ_W
  N_EventosCEvNS_NSIXe/    pipeline RED-100 xenón
  N_EventosCEvNS_NSIAr/    pipeline RED-100 argón (fondo de ³⁹Ar)
  chi2_nsi_generic/        motor NSI/sin²θ_W genérico (Xe/Ar/CONUS+)
python/                    figuras, tablas, validaciones, tablas NEST
doc/                       doc/metodologia.tex y figuras del pipeline
datos/                     entradas y salidas del pipeline (la mayoría se
                           regenera; ver .gitignore)
```

Cada carpeta de `FORTRAN90/` se compila por separado; no hay Makefile.
Detalle de la rama de argón en
`FORTRAN90/N_EventosCEvNS_NSIAr/README_Ar.txt`.

## Correr todo

```bash
./correr_todo.sh          # compila, corre el pipeline, genera figuras y tablas
./correr_todo.sh figuras  # solo regenera figuras/tablas
```

Dos cosas a saber:

- `nest.py`/`nest_Ar.py` usan semilla fija (`nestpy.RandomGen`), así que el
  resultado es reproducible bit a bit.
- El script no regenera `datos/chi2_sin2theta.dat`: `chi2.f90` (no se toca)
  tiene un bug de variable sin inicializar que a veces lo rompe. Queda
  congelado, ya validado.

## Compilar a mano

Orden de módulos obligatorio. Ejemplo, xenón:

```bash
cd FORTRAN90/N_EventosCEvNS_NSIXe
MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal        $MODS chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o chi2red100_nest   $MODS chi2red100_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest   $MODS mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest       $MODS red100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100PE          $MODS mod_detector.f90 red100PE.f90
gfortran -O2 -ffree-line-length-none -o chi2              constants.f90 chi2.f90
```

Argón: mismo patrón, sin `mod_detector.f90`/`red100PE.f90` (el argón no
tiene rama en fotoelectrones):

```bash
cd FORTRAN90/N_EventosCEvNS_NSIAr
MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal      $MODS chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o chi2_bkg        $MODS chi2_bkg_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest $MODS mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest     $MODS red100_nest.f90
```

`chi2_nsi_generic/`, `N_EventosCEvNS/` y `N_EventosCEvNS_NSI/` se compilan
igual, un binario por programa principal, módulos compartidos antes del
programa. Binarios/`.mod`/`.o` no se versionan (ver `.gitignore`).

## Paso a paso (lo que hace `correr_todo.sh`)

1. Compilar (sección anterior).
2. `python python/nest.py` y `python python/nest_Ar.py` (tablas NEST).
3. Correr los binarios (cada uno documenta su entrada/salida en su cabecera).
4. Correr los scripts de `python/` que grafican (`MPLBACKEND=Agg`). Estilo
   visual único en `python/estilo_tesis.py`.

El apéndice de `doc/metodologia.tex` tiene la tabla completa
programa → entrada → salida.

## Convenciones

- Todo parámetro se rastrea a una fuente publicada o a la salida del código.
- El fondo de ³⁹Ar (Ar) es una constante fija, no un parámetro libre.
- La incertidumbre se explora con corridas de reemplazo discretas
  (`python/sensitivity_scan.py`), no con un nuisance continuo.
