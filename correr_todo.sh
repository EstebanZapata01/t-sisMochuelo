#!/usr/bin/env bash
# =============================================================================
# correr_todo.sh -- reproduce TODO el pipeline de la tesis de un tirón:
#   1) tablas NEST (Python)          2) compila los binarios Fortran
#   3) corre los binarios Fortran    4) motor NSI/sin2theta genérico
#   5) genera las figuras y tablas (Python) citadas en doc/metodologia.tex
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
gfortran -O2 -ffree-line-length-none -o chi2red100_nest $MI chi2red100_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest $MI mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest     $MI red100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100PE        $MI mod_detector.f90 red100PE.f90
gfortran -O2 -ffree-line-length-none -o chi2            constants.f90 chi2.f90
cd "$ROOT"

echo "-- N_EventosCEvNS_NSIAr (Ar) --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSIAr"
MI="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
gfortran -O2 -ffree-line-length-none -o chi2_ideal      $MI chi2_ideal_nest.f90
gfortran -O2 -ffree-line-length-none -o chi2_bkg        $MI chi2_bkg_nest.f90
gfortran -O2 -ffree-line-length-none -o mainred100_nest $MI mainred100_nest.f90
gfortran -O2 -ffree-line-length-none -o red100_nest     $MI red100_nest.f90
cd "$ROOT"

echo "-- chi2_nsi_generic (motor NSI/sin2theta genérico) --"
cd "$ROOT/FORTRAN90/chi2_nsi_generic"
gfortran -O2 -ffree-line-length-none -o chi2_nsi_generic chi2_nsi_generic.f90
cd "$ROOT"

echo "-- N_EventosCEvNS_NSI (Ge/CONUS+ con NSI) --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSI"
MC="constants.f90 quenching.f90 xsections.f90 flux.f90 resolution.f90"
gfortran -O2 -ffree-line-length-none -o chi2_nsi_2D     $MC 2pchi2.f90
gfortran -O2 -ffree-line-length-none -o chi2_ON_OFF_1D  $MC chi2.f90
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
./chi2red100_nest
./mainred100_nest
./red100_nest
cd "$ROOT"

echo "-- Ar --"
cd "$ROOT/FORTRAN90/N_EventosCEvNS_NSIAr"
./chi2_ideal
./chi2_bkg
./mainred100_nest
./red100_nest
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
etapa "5/5 Figuras y tablas (Python, estilo_tesis.py)"
cd "$ROOT/python"
for s in ar_fondo_validacion blind_ring_XeAr chi2_perfil_generic chi2RED100 \
         helm_form_factor kappa_eee kappa_eee_plot ne_prediction_XeAr \
         recoil_spectrum_XeAr red100Xe sensibilidad_exposicion_completa \
         sensitivity_scan sin2theta_plot threshold_table expo_to_A90 \
         waterfall_XeAr valida_conus convergence_scan discovery_Z \
         tmax_kinematica dsigma_bare_XeAr nest_qy_fano_XeAr; do
    echo "-- $s.py --"
    $PY "$s.py"
done
cd "$ROOT"

# ----------------------------------------------------------------- resumen
etapa "Verificación: figuras/tablas citadas en doc/metodologia.tex"
FALTAN=0
grep -oE '\\includegraphics(\[[^]]*\])?\{[^}]+\}' "$ROOT/doc/metodologia.tex" \
    | grep -oE '\{[^}]+\}' | tr -d '{}' | sort -u | while read -r fig; do
  for ext in png pdf; do
    [ -f "$DATOS/$fig" ] && continue 2
    [ -f "$DATOS/$fig.$ext" ] && continue 2
  done
  echo "  FALTA: $fig"
  FALTAN=1
done
grep -oE 'tabla_[a-zA-Z0-9_]+\.tex' "$ROOT/doc/metodologia.tex" | sort -u | while read -r tab; do
  [ -f "$DATOS/$tab" ] || echo "  FALTA: $tab"
done

echo
echo "Listo. Figuras y tablas en $DATOS/."
