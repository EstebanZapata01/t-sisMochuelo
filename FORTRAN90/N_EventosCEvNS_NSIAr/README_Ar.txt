Programas del argon (misma estructura que N_EventosCEvNS_NSIXe)

La tabla de LArNEST se genera con python/nest_Ar.py.

Compilar desde esta carpeta (los modulos van en este orden):

  MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
  gfortran -O2 -o chi2_ideal      $MODS chi2_ideal_nest.f90
  gfortran -O2 -o chi2red100      $MODS chi2red100_nest.f90
  gfortran -O2 -o chi2_bkg        $MODS chi2_bkg_nest.f90
  gfortran -O2 -o mainred100_ar   $MODS mainred100_nest.f90
  gfortran -O2 -o red100_ar       $MODS red100_nest.f90
  gfortran -O2 -o cierre_pipeline $MODS cierre_pipeline.f90
  gfortran -O2 -o sigma_total     constants.f90 xsections_nest.f90 sigma_total.f90

Las salidas se escriben en ../../datos/.
