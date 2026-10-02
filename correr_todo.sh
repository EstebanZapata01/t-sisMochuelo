#!/usr/bin/env bash
# =============================================================================
# correr_todo.sh -- reproduce el pipeline de principio a fin:
#   1) tablas NEST (Python)          2) compila los binarios Fortran
#   3) corre los binarios Fortran    4) motor NSI/sin2theta generico
#   5) validaciones (Python) contra los resultados publicados de RED-100
#
# Uso:  ./correr_todo.sh            (todo)
#       ./correr_todo.sh validacion (solo la etapa 5, asume que 1-4 ya corrieron)
# =============================================================================
set -e
cd "$(dirname "$0")"
ROOT="$(pwd)"
DATOS="$ROOT/datos"
PY="python3"
export MPLBACKEND=Agg
mkdir -p "$DATOS"

etapa() { echo; echo "=== $1 ==="; }

SOLO_VALIDACION=0
[ "$1" = "validacion" ] && SOLO_VALIDACION=1

# chi2.f90 (Xe) y las validaciones contra RED-100 (valida_tabla1.py,
# chi2_perfil_generic.py, fig_perfiles_xe.py) necesitan el residuo ON-OFF y
# el perfil chi2(A) reales, digitalizados del PDF del paper (arXiv:2411.18641)
# con python/extrae_pdf_red100.py; se incluyen ya digitalizados en datos/.
for f in red100_residuo_ONOFF_fig8.csv red100_perfil_chi2_fig9.csv; do
    [ -f "$DATOS/$f" ] || { echo "ERROR: falta datos/$f (deberia venir con el repo)."; exit 1; }
done

# ----------------------------------------------------------------- etapa 1
if [ "$SOLO_VALIDACION" = 0 ]; then
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

echo "-- chi2_nsi_generic (motor NSI/sin2theta generico) --"
cd "$ROOT/FORTRAN90/chi2_nsi_generic"
gfortran -O2 -ffree-line-length-none -o chi2_nsi_generic chi2_nsi_generic.f90
cd "$ROOT"

echo "-- N_EventosCEvNS_NSI (Ge/CONUS+ con NSI) --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSI"
MC="constants.f90 quenching.f90 xsections.f90 flux.f90 resolution.f90"
gfortran -O2 -ffree-line-length-none -o chi2_nsi_2D     $MC 2pchi2.f90
gfortran -O2 -ffree-line-length-none -o eventos_conus   $MC main.f90
cd "$ROOT"

echo "-- N_EventosCEvNS (Ge/CONUS+ original, validacion sin2theta) --"
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
./chi2              # lee ionization_spectra_detallado.dat (de red100PE); escribe generic_input_Xe.dat
./mainred100_nest
./red100_nest
./cierre_pipeline
./sigma_total
USE_HELM=1 ./sigma_total
cd "$ROOT"

echo "-- Ar --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSIAr"
./chi2_ideal
./chi2_bkg
./mainred100_nest
./red100_nest
./cierre_pipeline
./sigma_total
USE_HELM=1 ./sigma_total
cd "$ROOT"

echo "-- Ge/CONUS+ --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS"
./eventos_conus
cd "$ROOT"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSI"
./chi2_nsi_2D           # escribe generic_input_conus.dat (validacion cruzada)
cd "$ROOT"

# ----------------------------------------------------------------- etapa 4
etapa "4/5 Motor NSI/sin2theta generico (validacion cruzada Xe/CONUS+)"
cd "$ROOT/FORTRAN90/chi2_nsi_generic"
./chi2_nsi_generic "$DATOS/generic_input_Xe.dat" "$DATOS/generic_Xe"
./chi2_nsi_generic "$DATOS/generic_input_conus.dat" "$DATOS/generic_conus"
cd "$ROOT"
fi

# ----------------------------------------------------------------- etapa 5
etapa "5/5 Validaciones (Python) contra los resultados publicados de RED-100"
cd "$ROOT/python"
for s in valida_tabla1 chi2_perfil_generic fig_perfiles_xe roi_cuts_validacion_Xe; do
    echo "-- $s.py --"
    $PY "$s.py"
done
cd "$ROOT"

echo
echo "Listo. Resultados en $DATOS/."
