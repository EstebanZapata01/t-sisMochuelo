!=======================================================================
! Programa: main_red100 (Versión Argón - ROI Extendida)
!   Genera el Asimov Dataset (Eventos totales esperados SM)
!   para RED-100 (Ar) en la ROI discreta extendida (1 a 20 electrones extraídos).
!=======================================================================
program main_red100
  use constants
  use quenching
  use xsections
  use flux, only: kopeikin_spectrum, mueller_spectrum, fission_frac, &
                  E_nu_min_flux, E_nu_max
  implicit none

  ! --- Definición de la ROI ---
  integer, parameter :: bin_min = 1
  integer, parameter :: bin_max = 20

  integer, parameter :: n_T = 500
  integer, parameter :: n_E = 2000
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, integrando
  real(dp) :: tasa_Comb, QW_SM
  integer :: i_T, i_E, n_bin, u_out
  
  ! Arreglo ajustado a la ROI dinámica
  real(dp) :: R_SM(bin_min:bin_max)
  
  real(dp) :: E_ee, n_creados, n_extraidos

  ! Parámetros de Exposición y Reactor (Kalinin)
  real(dp), parameter :: P_th = 3.1_dp
  real(dp), parameter :: E_f_MeV = 204.0_dp
  real(dp), parameter :: L_m = 19.0_dp
  real(dp), parameter :: dias_exposicion = 331.0_dp
  real(dp) :: fission_rate, flux_factor, L_cm, E_f_joules
  real(dp) :: atoms_per_kg, sec_per_day, conv_keV

  character(len=200) :: outdir, filename

  ! Ajusta esta ruta a tu carpeta
  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  filename = trim(outdir) // 'asimov_ar.dat'

  ! --- Carga débil SM ---
  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*0.23857_dp)/2.0_dp * Z_Ge
  QV2 = QW_SM**2

  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp
  conv_keV     = 1000.0_dp

  ! --- Factor de flujo macroscópico ---
  E_f_joules   = E_f_MeV * 1.60218d-13
  L_cm         = L_m * 100.0_dp
  fission_rate = (P_th * 1.0d9) / E_f_joules
  flux_factor  = fission_rate / (4.0_dp * PI * L_cm**2)

  ! Inicializar el array de la ROI en ceros
  R_SM = 0.0_dp

  ! OJO: Para llegar a 20 electrones en Argón (con eta=23.6eV), 
  ! necesitas retrocesos de casi 10 keV. Por seguridad, aumentamos T_nr_max
  T_nr_min = 0.1_dp / 1000.0_dp
  T_nr_max = 12.0_dp / 1000.0_dp  ! Ampliado de 3 a 12 keVnr
  dT = (T_nr_max - T_nr_min) / (n_T - 1)
  dT_keV = dT * 1000.0_dp
  dE = E_nu_max / (n_E - 1)

  write(*,*) '============================================'
  write(*,*) '    GENERADOR ASIMOV RED-100 (Ar) S2-Only   '
  write(*,*) '============================================'
  write(*,*) 'Masa Activa (kg): ', total_mass_kg
  write(*,*) 'Exposición (días): ', dias_exposicion
  write(*,*) 'Eficiencia Extracción (EEE): ', EEE
  write(*,*) 'Rango de Integración (keVnr): ', T_nr_min*1000.0_dp, ' a ', T_nr_max*1000.0_dp
  write(*,*) '============================================'

  ! Bucle principal termodinámico
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Comb = 0.0_dp

     ! Integración sobre el flujo de neutrinos
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle

        if (E_nu < E_nu_min_flux) then
           integrando = kopeikin_spectrum(E_nu) * dsigma_dT(E_nu, T_nr)
        else
           integrando = ( fission_frac(1)*mueller_spectrum(E_nu,1) &
                        + fission_frac(2)*mueller_spectrum(E_nu,2) &
                        + fission_frac(3)*mueller_spectrum(E_nu,3) &
                        + fission_frac(4)*mueller_spectrum(E_nu,4) ) &
                        * dsigma_dT(E_nu, T_nr)
        end if

        if (i_E == 1 .or. i_E == n_E) then
           tasa_Comb = tasa_Comb + 0.5_dp * integrando * dE
        else
           tasa_Comb = tasa_Comb + integrando * dE
        end if
     end do

     ! Tasa diferencial en eventos / (keV * kg * día)
     tasa_Comb = tasa_Comb * flux_factor * conv * atoms_per_kg * sec_per_day / conv_keV

     ! Transformación de retroceso a conteo discreto de electrones
     E_ee = T_nr * QF(T_nr)
     n_creados = E_ee / eta
     n_extraidos = n_creados * EEE
     n_bin = nint(n_extraidos)

     ! Acumular EVENTOS TOTALES dinámicamente
     if (n_bin >= bin_min .and. n_bin <= bin_max) then
        R_SM(n_bin) = R_SM(n_bin) + (tasa_Comb * dT_keV * total_mass_kg * dias_exposicion)
     end if
  end do

  ! Imprimir y guardar resultados
  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# N_e   Eventos_Totales'
  do n_bin = bin_min, bin_max
     write(u_out, '(I2, 2X, ES15.6)') n_bin, R_SM(n_bin)
     write(*, '(A, I2, A, F12.4)') 'Bin Ne = ', n_bin, ' -> Eventos totales = ', R_SM(n_bin)
  end do
  close(u_out)

  write(*,*) '============================================'
  write(*,*) 'Archivo guardado: ', trim(filename)

end program main_red100
