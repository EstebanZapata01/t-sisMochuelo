#!/usr/bin/env bash
# =============================================================================
# correr_todo.sh -- reproduce TODO el pipeline de la tesis de un tirón:
#   1) tablas NEST (Python)          2) compila los binarios Fortran
#   3) corre los binarios Fortran    4) motor NSI/sin2theta genérico
#   5) regenera las figuras de la tesis (python/figuras_tesis.txt) y las
#      copia a Proyecto_LaTeX/figuras_metodologia
#
# Uso:  ./correr_todo.sh          (todo)
#       ./correr_todo.sh figuras  (solo la etapa 5, asume que 1-4 ya corrieron)
#
# El detalle de cada comando de compilación está documentado también en
# README.md; este script simplemente los encadena en el orden correcto.
# =============================================================================
set -e
cd "$(dirname "$0")"
ROOT="$(pwd)"
DATOS="$ROOT/datos"
PY="python3"
export MPLBACKEND=Agg

etapa() { echo; echo "=== $1 ==="; }

SOLO_FIGURAS=0
[ "$1" = "figuras" ] && SOLO_FIGURAS=1

# ----------------------------------------------------------------- etapa 1
if [ "$SOLO_FIGURAS" = 0 ]; then
etapa "1/5 Tablas NEST (insumo de Fortran)"
cd "$ROOT/python"
$PY nest.py
$PY nest_Ar.py
cd "$ROOT"

# ----------------------------------------------------------------- etapa 2
etapa "2/5 Compilando binarios Fortran"

echo "-- N_EventosCEvNS_NSIXe (Xe, RED-100) --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSIXe"
MI="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal      $MI chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest $MI mod_detector.f90 mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest     $MI red100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100PE        $MI mod_detector.f90 red100PE.f90
gfortran -O2 -ffree-line-length-none -o chi2            constants.f90 mod_detector.f90 chi2.f90
gfortran -O2 -ffree-line-length-none -o cierre_pipeline $MI mod_detector.f90 cierre_pipeline.f90
gfortran -O2 -ffree-line-length-none -o sigma_total     $MI sigma_total.f90
cd "$ROOT"

echo "-- N_EventosCEvNS_NSIAr (Ar) --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSIAr"
MI="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal      $MI chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o chi2_bkg        $MI chi2_bkg_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest $MI mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest     $MI red100_nest.f90
gfortran -O2 -ffree-line-length-none -o cierre_pipeline $MI cierre_pipeline.f90
gfortran -O2 -ffree-line-length-none -o sigma_total     $MI sigma_total.f90
cd "$ROOT"

echo "-- chi2_nsi_generic (motor NSI/sin2theta genérico) --"
cd "$ROOT/FORTRAN90/chi2_nsi_generic"
gfortran -O2 -ffree-line-length-none -o chi2_nsi_generic chi2_nsi_generic.f90
cd "$ROOT"

echo "-- N_EventosCEvNS_NSI (Ge/CONUS+ con NSI) --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSI"
MC="constants.f90 quenching.f90 xsections.f90 flux.f90 resolution.f90"
gfortran -O2 -ffree-line-length-none -o chi2_nsi_2D     $MC 2pchi2.f90
gfortran -O2 -ffree-line-length-none -o eventos_conus   $MC main.f90
cd "$ROOT"

echo "-- N_EventosCEvNS (Ge/CONUS+ original, validación sin2theta) --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS"
MC="constants.f90 quenching.f90 xsections.f90 flux.f90 resolution.f90"
gfortran -O2 -ffree-line-length-none -o chi2_nsi_2D     $MC 2pchi2.f90
gfortran -O2 -ffree-line-length-none -o chi2_sin2theta  $MC chi2.f90
gfortran -O2 -ffree-line-length-none -o eventos_conus   $MC main.f90
cd "$ROOT"

# ----------------------------------------------------------------- etapa 3
etapa "3/5 Corriendo los binarios (genera los .dat base)"

echo "-- Xe --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSIXe"
./chi2_ideal
./red100PE
./chi2                 # lee ionization_spectra_detallado.dat (de red100PE); escribe generic_input_Xe.dat
./mainred100_nest
./red100_nest
./cierre_pipeline      # datos/cierre_{binomial,dsigma,factorizacion}_Xe.dat (cierres, rutinas reales)
./sigma_total          # datos/sigma_total_Xe.dat (Helm, default)
USE_HELM=0 ./sigma_total   # datos/sigma_total_bare_Xe.dat (F^2=1)
cd "$ROOT"

echo "-- Ar --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSIAr"
./chi2_ideal
./chi2_bkg
./mainred100_nest
./red100_nest
./cierre_pipeline      # datos/cierre_{binomial,dsigma,factorizacion}_Ar.dat
./sigma_total          # datos/sigma_total_Ar.dat (Helm, default)
USE_HELM=0 ./sigma_total   # datos/sigma_total_bare_Ar.dat (F^2=1)
cd "$ROOT"

echo "-- Ge/CONUS+ --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS"
# NO se corre ./chi2_sin2theta aqui: bug conocido en chi2.f90, ver README.md.
./eventos_conus
cd "$ROOT"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSI"
./chi2_nsi_2D           # escribe generic_input_conus.dat (validación cruzada)
cd "$ROOT"

# ----------------------------------------------------------------- etapa 4
etapa "4/5 Motor NSI/sin2theta genérico (validación cruzada Xe/CONUS+)"
cd "$ROOT/FORTRAN90/chi2_nsi_generic"
./chi2_nsi_generic "$DATOS/generic_input_Xe.dat"    "$DATOS/generic_Xe"
./chi2_nsi_generic "$DATOS/generic_input_conus.dat" "$DATOS/generic_conus"
cd "$ROOT"
fi

# ----------------------------------------------------------------- etapa 5
etapa "5/5 Figuras de la tesis (python/figuras_tesis.txt) -> Proyecto_LaTeX"
cd "$ROOT/python"
TESIS="$ROOT/../Proyecto_LaTeX/figuras_metodologia"
while read -r script figura; do
    [ -z "$script" ] && continue
    case "$script" in \#*) continue ;; esac
    echo "-- $script.py --"
    $PY "$script.py"
    cp "$DATOS/$figura" "$TESIS/"
done < figuras_tesis.txt
cd "$ROOT"

echo
echo "Listo. Figuras en $DATOS/ y en $ROOT/../Proyecto_LaTeX/figuras_metodologia/."
