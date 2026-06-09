!=======================================================================
! Programa: red100 Ar (Versión Optimizada Tesis)
!   Cálculo fenomenológico riguroso para el detector RED-100 (Argón)
!   1) Espectro de retroceso nuclear (dN/dT_nr) con normalización absoluta
!   2) Espectro en electrones de ionización conservando probabilidad total
!=======================================================================
program red100
  use constants
  use quenching
  use xsections
  use flux, only: kopeikin_spectrum, mueller_spectrum, &
                  fission_frac, E_nu_min_flux, E_nu_max
  implicit none

  ! Parámetros de malla topológica
  integer, parameter :: n_T = 500        ! Nodos en T_nr
  integer, parameter :: n_E = 2000       ! Nodos en E_nu
  integer, parameter :: n_ion = 20       ! Límite de electrones observables

  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, QW_SM
  real(dp) :: integrando_kop, integrando_mue, integrando_comb
  integer :: i_T, i_E, i_bin, u_out
  real(dp) :: tasa_Kop, tasa_Mue, tasa_Comb
  character(len=200) :: outdir, filename

  ! Parámetros físicos del Reactor Kalinin
  real(dp), parameter :: P_th = 3.1_dp
  real(dp), parameter :: E_f_MeV = 204.0_dp
  real(dp), parameter :: L_m = 19.0_dp
  real(dp) :: fission_rate, flux_factor, L_cm, E_f_joules

  ! Factores de conversión volumétrica y temporal
  real(dp) :: atoms_per_kg, sec_per_day, conv_keV

  ! Estructuras de datos para el mapeo a ionización
  real(dp), allocatable :: array_tasa_Comb(:)
  real(dp), allocatable :: tasa_ion_creados(:), tasa_ion_extraidos(:)
  real(dp) :: E_ee, n_creados_medio, n_extraidos_medio, peso
  integer :: n_bin_c, n_bin_e

  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'

  ! --- Carga débil SM y constantes estequiométricas ---
  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*0.23857_dp)/2.0_dp * Z_Ge
  QV2 = QW_SM**2

  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp
  conv_keV     = 1000.0_dp

  ! --- Cálculo del factor de flujo macroscópico ---
  E_f_joules   = E_f_MeV * 1.60218d-13
  L_cm         = L_m * 100.0_dp
  fission_rate = (P_th * 1.0d9) / E_f_joules
  flux_factor  = fission_rate / (4.0_dp * PI * L_cm**2)
  
  ! --- Asignación de memoria dinámica ---
  allocate(array_tasa_Comb(n_T))
  allocate(tasa_ion_creados(n_ion), tasa_ion_extraidos(n_ion))
  tasa_ion_creados(:) = 0.0_dp
  tasa_ion_extraidos(:) = 0.0_dp

  write(*,*) '=== DIAGNÓSTICO red100 (v8.3 FINAL ARGÓN) ==='
  write(*,*) 'Flux factor (cm^-2 s^-1) = ', flux_factor
  write(*,*) '============================================='

  ! ==================== 1. Espectros de retroceso (Espacio Continuo) ====================
  T_nr_min = 0.1_dp / 1000.0_dp
  T_nr_max = 12.0_dp / 1000.0_dp  ! 12 keV para abarcar toda la cinemática del Argón
  dT = (T_nr_max - T_nr_min) / (n_T - 1)
  dE = E_nu_max / (n_E - 1)

  filename = trim(outdir) // 'recoil_spectraAr.dat'
  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# T_nr(keV)   Kop_full   Mue_pure   Combinado'

  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Kop  = 0.0_dp
     tasa_Mue  = 0.0_dp
     tasa_Comb = 0.0_dp

     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle

        ! 1. Calcular SIEMPRE Kopeikin para obtener su curva completa
        integrando_kop = kopeikin_spectrum(E_nu) * dsigma_dT(E_nu, T_nr)

        ! 2. Calcular SIEMPRE Mueller para obtener su curva completa
        integrando_mue = ( fission_frac(1)*mueller_spectrum(E_nu,1) &
                         + fission_frac(2)*mueller_spectrum(E_nu,2) &
                         + fission_frac(3)*mueller_spectrum(E_nu,3) &
                         + fission_frac(4)*mueller_spectrum(E_nu,4) ) &
                         * dsigma_dT(E_nu, T_nr)

        ! 3. Construir el flujo COMBINADO oficial (Kopeikin < 2 MeV, Mueller >= 2 MeV)
        if (E_nu < E_nu_min_flux) then
           integrando_comb = integrando_kop
        else
           integrando_comb = integrando_mue
        end if

        ! Aplicar regla del trapecio para la integración en energía del neutrino
        if (i_E == 1 .or. i_E == n_E) then
           peso = 0.5_dp
        else
           peso = 1.0_dp
        end if

        tasa_Kop  = tasa_Kop  + peso * integrando_kop * dE
        tasa_Mue  = tasa_Mue  + peso * integrando_mue * dE
        tasa_Comb = tasa_Comb + peso * integrando_comb * dE
     end do

     ! Escalado físico final a eventos / (keV kg día)
     tasa_Kop  = tasa_Kop  * flux_factor * conv * atoms_per_kg * sec_per_day / conv_keV
     tasa_Mue  = tasa_Mue  * flux_factor * conv * atoms_per_kg * sec_per_day / conv_keV
     tasa_Comb = tasa_Comb * flux_factor * conv * atoms_per_kg * sec_per_day / conv_keV

     ! Almacenamiento en memoria para la transformación discreta posterior
     array_tasa_Comb(i_T) = tasa_Comb

     write(u_out, '(F10.4, 3ES15.6)') T_nr*1.0e3_dp, tasa_Kop, tasa_Mue, tasa_Comb
  end do
  close(u_out)
  write(*,*) 'Archivo recoil_spectraAr.dat generado exitosamente.'

  ! ==================== 2. Espectro de ionización (Espacio Discreto) ====================
  filename = trim(outdir) // 'ionization_spectraAr.dat'
  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# N_e   creados   extraídos'

  ! Cuadratura estricta con factor de Lindhard
  dT_keV = dT * 1000.0_dp
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     
     ! Aplicación del Quenching Factor (QF) analítico
     E_ee = T_nr * QF(T_nr) 
     
     ! Determinación de valores esperados de electrones
     n_creados_medio = E_ee / eta
     n_extraidos_medio = n_creados_medio * EEE
     
     ! Mapeo estocástico mediante redondeo al entero más cercano
     n_bin_c = nint(n_creados_medio)
     n_bin_e = nint(n_extraidos_medio)
     
     ! Aplicar peso trapezoidal riguroso para la integración en T_nr
     if (i_T == 1 .or. i_T == n_T) then
         peso = 0.5_dp
     else
         peso = 1.0_dp
     end if
     
     ! Acumulación discreta conservando la probabilidad total (integrando el Combinado)
     if (n_bin_c >= 1 .and. n_bin_c <= n_ion) then
        tasa_ion_creados(n_bin_c) = tasa_ion_creados(n_bin_c) + array_tasa_Comb(i_T) * dT_keV * peso
     end if
     
     if (n_bin_e >= 1 .and. n_bin_e <= n_ion) then
        tasa_ion_extraidos(n_bin_e) = tasa_ion_extraidos(n_bin_e) + array_tasa_Comb(i_T) * dT_keV * peso
     end if
  end do

  ! Exportación de la distribución discretizada final
  do i_bin = 1, n_ion
     write(u_out, '(F6.1, 2ES15.6)') real(i_bin, dp), tasa_ion_creados(i_bin), tasa_ion_extraidos(i_bin)
  end do
  close(u_out)
  write(*,*) 'Archivo ionization_spectraAr.dat generado.'

  ! Liberación de memoria
  deallocate(array_tasa_Comb, tasa_ion_creados, tasa_ion_extraidos)
  write(*,*) '=== COMPILACIÓN Y EJECUCIÓN FINALIZADA ==='

end program red100
