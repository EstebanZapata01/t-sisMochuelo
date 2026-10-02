# RED-100: sensibilidad a NSI y a sin²θ_W, Xe vs. Ar

Motor de calculo de la tesis de grado de Esteban Zapata (Universidad de
Pamplona), que compara la sensibilidad de un detector estilo RED-100
([arXiv:2411.18641](https://arxiv.org/abs/2411.18641), CEνNS de
antineutrinos de reactor) a interacciones no estandar (NSI) y al angulo de
mezcla debil sin²θ_W, para un blanco de xenon (validado contra el dato real
2024 de la colaboracion) y uno de argon (proyeccion, con los parametros
publicados en D. Akimov et al., Physics 5, 492 (2023) y en
[arXiv:2510.16404](https://arxiv.org/abs/2510.16404)).

Este repositorio tiene el motor de fisica (Fortran) y solo las validaciones
que lo contrastan directamente contra los resultados publicados de RED-100
(Python). No incluye los scripts que generan las demas figuras de la tesis,
ni los analisis de sensibilidad/discusion, ni el documento de la tesis en
si: solo el calculo y su verificacion contra el paper.

## Estructura

```
FORTRAN90/
  N_EventosCEvNS/          Ge/CONUS+ original (sin NSI), referencia de validacion
  N_EventosCEvNS_NSI/      Ge/CONUS+ con NSI 2D y barrido de sin²θ_W
  N_EventosCEvNS_NSIXe/    pipeline RED-100, blanco de xenon
  N_EventosCEvNS_NSIAr/    pipeline RED-100, blanco de argon
  chi2_nsi_generic/        motor NSI/sin²θ_W generico (el mismo codigo
                            corre con los datos de Xe y de CONUS+)
python/
  nest.py, nest_Ar.py        tablas NEST/LArNEST (insumo de Fortran)
  extrae_pdf_red100.py       digitaliza del PDF del paper el residuo ON-OFF
                              y el perfil chi2(A) (ya incluidos en datos/,
                              ver "Correr todo")
  valida_tabla1.py           reproduce la Tabla I de RED-100 (SM2018) con
                              SOLO datos publicados por la colaboracion
  roi_cuts_validacion_Xe.py  valida la reconstruccion en PE/ROI contra la
                              Fig. 6 de RED-100 (senal antes/despues de cortes)
  chi2_perfil_generic.py     prueba visual de que el mismo binario (Xe y
                              CONUS+) reproduce ambos perfiles publicados
  fig_perfiles_xe.py         perfil Δχ²(A) de xenon: motor + señal de la
                              colaboracion vs. cadena + predicción propia
  nsi_core.py, leer_fortran.py, estilo_tesis.py
                              modulos de soporte (no se corren solos)
datos/                     entradas/salidas del pipeline; solo el residuo
                            ON-OFF y el perfil chi2(A) digitalizados del
                            paper estan versionados (no regenerables sin el
                            PDF); el resto se regenera (ver .gitignore)
```

Cada carpeta de `FORTRAN90/` se compila por separado; no hay Makefile.
`FORTRAN90/N_EventosCEvNS_NSIAr/README_Ar.txt` detalla que cambia entre la
rama de Xe y la de Ar.

## Correr todo

```bash
./correr_todo.sh             # compila, corre el pipeline, corre las validaciones
./correr_todo.sh validacion  # solo la etapa 5 (Python), asume 1-4 ya corridas
```

Dos cosas a saber:

- `nest.py`/`nest_Ar.py` usan semilla fija (`nestpy.RandomGen`), asi que el
  resultado es reproducible bit a bit. Necesitan `nestpy` instalado.
- `chi2.f90` (el ajuste real de xenon contra el residuo ON-OFF publicado) y
  las validaciones que dependen de el (`valida_tabla1.py`,
  `chi2_perfil_generic.py`, `fig_perfiles_xe.py`) necesitan el residuo
  ON-OFF y el perfil chi2(A) reales de RED-100, digitalizados del PDF del
  paper con `python/extrae_pdf_red100.py`. Esas dos tablas YA estan
  incluidas (`datos/red100_residuo_ONOFF_fig8.csv`,
  `datos/red100_perfil_chi2_fig9.csv`: no se pueden regenerar sin el PDF,
  no incluido aqui); `correr_todo.sh` falla rapido si no estan.

## Compilar a mano

Orden de modulos obligatorio. Ejemplo, xenon:

```bash
cd FORTRAN90/N_EventosCEvNS_NSIXe
MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal        $MODS chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest   $MODS mod_detector.f90 mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest       $MODS red100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100PE          $MODS mod_detector.f90 red100PE.f90
gfortran -O2 -ffree-line-length-none -o chi2              constants.f90 mod_detector.f90 chi2.f90
```

Argon: mismo patron, sin `mod_detector.f90`/`red100PE.f90` (el argon no
tiene rama en fotoelectrones, ver `README_Ar.txt`):

```bash
cd FORTRAN90/N_EventosCEvNS_NSIAr
MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal      $MODS chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o chi2_bkg        $MODS chi2_bkg_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest $MODS mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest     $MODS red100_nest.f90
```

`chi2_nsi_generic/`, `N_EventosCEvNS/` y `N_EventosCEvNS_NSI/` se compilan
igual: un binario por programa principal, modulos compartidos antes del
programa. Binarios, `.mod` y `.o` no se versionan (ver `.gitignore`).

## Convenciones

- Todo parametro numerico se rastrea a una fuente publicada (citada en el
  comentario de cabecera del archivo que lo define) o a la salida de otra
  etapa del pipeline.
- El fondo de Ar es una constante fija, el nivel declarado por la
  referencia publicada (S/√B≈4 a 62 kg·dia, plano en N_e=1..5); no se
  simula ningun proceso de fondo ni se trata como parametro libre.
- Las validaciones comparan contra datos y cifras publicadas sin
  modificarlos; cuando un resultado propio difiere del publicado, el
  comentario de cabecera del programa correspondiente explica por que.
