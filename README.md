# RED-100: sensibilidad a NSI y a sin²θ_W, Xe vs. Ar

Pipeline de la tesis de grado de Esteban Zapata (Universidad de Pamplona).
Compara la sensibilidad de un detector estilo RED-100 (arXiv:2411.18641,
CEνNS de antineutrinos de reactor) a NSI y al ángulo de mezcla débil, para
un blanco de xenón (validado contra el dato real 2024) y uno de argón
(proyección, ref. [46] + ReD 2025).

Repositorio PRIVADO (área de trabajo completa, con flujo visual). La
metodología completa vive en `Proyecto_LaTeX/` (tesis); el subconjunto que
sale a GitHub (solo la cadena de cálculo Fortran) se arma aparte, sin este
README ni las herramientas de la sección "Flujo de trabajo".

## Estructura

```
FORTRAN90/
  N_EventosCEvNS/          Ge/CONUS+ original (sin NSI)
  N_EventosCEvNS_NSI/      Ge/CONUS+ con NSI 2D y barrido de sin²θ_W
  N_EventosCEvNS_NSIXe/    pipeline RED-100 xenón
  N_EventosCEvNS_NSIAr/    pipeline RED-100 argón (fondo = nivel declarado ref.[46])
  chi2_nsi_generic/        motor NSI/sin²θ_W genérico (Xe/Ar/CONUS+)
python/
  nest.py, nest_Ar.py      tablas NEST (nestpy), insumo de Tnr_to_e.f90
  figuras_tesis.txt        lista script -> figura, en el orden de la tesis
  fig_*.py, otros          las figuras y cifras que citan la tesis y la
                           charla (ver cabecera de cada script: entrada/salida)
  estilo_tesis.py          estilo visual único; nsi_core.py, leer_fortran.py
                           son los otros módulos compartidos
  vigilar.py, explorar.ipynb   flujo de trabajo visual (ver abajo)
datos/, datos_tutor/       entradas y salidas del pipeline (la mayoría se
                           regenera; ver .gitignore)
```

Cada carpeta de `FORTRAN90/` se compila por separado; no hay Makefile.
Detalle de la rama de argón en
`FORTRAN90/N_EventosCEvNS_NSIAr/README_Ar.txt`.

## Correr todo

```bash
./correr_todo.sh          # compila, corre el pipeline, regenera las figuras de la tesis
./correr_todo.sh figuras  # solo la etapa 5 (figuras), asume que 1-4 ya corrieron
```

`./correr_todo.sh figuras` deja las 21 figuras de `python/figuras_tesis.txt`
en `datos/` y copiadas en `Proyecto_LaTeX/figuras_metodologia/`.

Dos cosas a saber:

- `nest.py`/`nest_Ar.py` usan semilla fija (`nestpy.RandomGen`), así que el
  resultado es reproducible bit a bit.
- El script no regenera `datos/chi2_sin2theta.dat`: `chi2.f90` (no se toca)
  tiene un bug de variable sin inicializar que a veces lo rompe. Queda
  congelado, ya validado.

## Flujo de trabajo (privado)

Para ver en segundos qué pasa con una figura al tocar un script o una
tabla, sin correr `correr_todo.sh` entero:

```bash
cd python && python3 vigilar.py
```

Mira las fechas de `python/*.py` y `datos/*.dat` cada segundo; si cambia un
script de `figuras_tesis.txt` (o un módulo que alguno de ellos importa, o
cualquier `.dat`), lo vuelve a correr y copia su PNG a
`Proyecto_LaTeX/figuras_metodologia/`. El log de cada corrida queda en
`datos/logs/<script>.txt`.

`python/explorar.ipynb` tiene tres paneles interactivos (ipywidgets) para
mover un parámetro (exposición, blanco, umbral, par NSI, t_ON/t_OFF) y ver
$A_{90}$ o el plano NSI al instante, reutilizando las mismas fórmulas de
`fig_sensibilidad_mapa.py`/`tabla_fondo_Ar.py`/`nsi_core.py` sin redefinir
física nueva.

## Compilar a mano

Orden de módulos obligatorio. Ejemplo, xenón:

```bash
cd FORTRAN90/N_EventosCEvNS_NSIXe
MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal        $MODS chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest   $MODS mod_detector.f90 mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest       $MODS red100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100PE          $MODS mod_detector.f90 red100PE.f90
gfortran -O2 -ffree-line-length-none -o chi2              constants.f90 mod_detector.f90 chi2.f90
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
4. Motor NSI/sin²θ_W genérico (validación cruzada Xe/CONUS+).
5. Correr los scripts de `python/figuras_tesis.txt` (`MPLBACKEND=Agg`) y
   copiar cada PNG a `Proyecto_LaTeX/figuras_metodologia/`. Estilo visual
   único en `python/estilo_tesis.py`.

## Convenciones

- Factor de forma de Helm (Lewin & Smith, Astropart. Phys. 6, 87 (1996))
  ACTIVADO por defecto en toda la cadena (Xe, Ar, CONUS+/Ge):
  `USE_HELM=0 ./binario` lo apaga (F²=1) para sensibilidad. `sigma_total`
  escribe `sigma_total_{Xe,Ar}.dat` con Helm y `sigma_total_bare_{Xe,Ar}.dat`
  con F²=1.
- Todo parámetro se rastrea a una fuente publicada o a la salida del código.
- El fondo de Ar es una constante fija (nivel declarado por la ref. [46],
  S/√B≈4 a 62 kg·día, plano en N_e=1..5), no un parámetro libre ni un
  fondo simulado.
