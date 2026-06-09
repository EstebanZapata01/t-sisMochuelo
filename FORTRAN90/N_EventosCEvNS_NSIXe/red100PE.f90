!=======================================================================
! Programa: red100PE_detallado (Plantillas PCHIP - Flujo Explícito)
!=======================================================================
program red100PE_detallado
  use constants
  use mod_tnr_to_e        
  use xsections_nest      
  use mod_stats           
  use mod_detector        
  use flux, only: flujo_diferencial, E_nu_max
  implicit none

  integer, parameter :: n_T = 500, n_E = 2000, n_ion = 15
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, peso
  integer :: i_T, i_E, i_bin, u_out, i_pe
  real(dp) :: tasa_Comb, S_pe, bin_width_pe, QW_SM
  real(dp), allocatable :: array_tasa_Comb(:), tasa_ion_extraidos(:), contribuciones(:)
  real(dp) :: total_bin, n_creados_medio, n_extraidos_medio
  character(len=250) :: outdir, filename, datafile
  real(dp) :: atoms_per_kg, sec_per_day, dummy

  bin_width_pe = 5.0_dp
  
  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  datafile = trim(outdir) // 'templates_SE_pchip.dat'
  filename = trim(outdir) // 'ionization_spectra_detallado.dat'

  call inicializar_nest()
  call cargar_SE_data(datafile)
  
  ! Disparamos la inicialización automática del flujo
  dummy = flujo_diferencial(1.0_dp)

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2

  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day = 86400.0_dp

  allocate(array_tasa_Comb(n_T), tasa_ion_extraidos(n_ion), contribuciones(6))
  tasa_ion_extraidos(:) = 0.0_dp

  T_nr_min = 0.1_dp / 1000.0_dp; T_nr_max = 3.0_dp / 1000.0_dp
  dT = (T_nr_max - T_nr_min) / (n_T - 1); dT_keV = dT * 1000.0_dp
  dE = E_nu_max / (n_E - 1)

  write(*,*) "-> Integrando el espectro continuo..."
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Comb = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle
        peso = merge(0.5_dp, 1.0_dp, i_E == 1 .or. i_E == n_E)
        tasa_Comb = tasa_Comb + (flujo_diferencial(E_nu) * dsigma_dT(E_nu, T_nr) * dE * peso)
     end do
     array_tasa_Comb(i_T) = tasa_Comb * atoms_per_kg * sec_per_day
  end do

  write(*,*) "-> Aplicando estadística de Poisson..."
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     n_creados_medio = obtener_electrones_creados(T_nr * 1000.0_dp)
     n_extraidos_medio = n_creados_medio * EEE
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)
     
     do i_bin = 1, n_ion
        tasa_ion_extraidos(i_bin) = tasa_ion_extraidos(i_bin) + &
     (array_tasa_Comb(i_T) * dT_keV * peso * livetime_frac * poisson_prob(i_bin, n_extraidos_medio))
     end do
  end do


  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# PE_center   Total   1SE   2SE   3SE   4SE   5SE   6SE'

  write(*,*) "-> Escribiendo formato PCHIP..."
  do i_pe = 1, 400
     S_pe = (real(i_pe, dp) - 0.5_dp) * bin_width_pe
     total_bin = 0.0_dp
     
     do i_bin = 1, 6
        contribuciones(i_bin) = tasa_ion_extraidos(i_bin) * respuesta_empirica(S_pe, i_bin) * bin_width_pe
        total_bin = total_bin + contribuciones(i_bin)
     end do
     
     write(u_out, '(F8.1, E15.6)', advance='no') S_pe, total_bin
     do i_bin = 1, 6
        write(u_out, '(E15.6)', advance='no') contribuciones(i_bin)
     end do
     write(u_out, *)
  end do
  close(u_out)

  write(*,*) "=== EXITOSO: Espectro guardado en Eventos/bin/kg/dia ==="
end program red100PE_detallado
