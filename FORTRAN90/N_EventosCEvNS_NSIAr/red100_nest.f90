!=======================================================================
! Programa: red100_nest (Poisson → Binomial, solo creados y extraídos)
!   CLON de la carpeta de Xe; rutas de salida *_Ar.dat.
! Pipeline: etapa "retroceso -> ionizacion" (diagnostico creados/extraidos).
! Tesis   : metodologia.tex Sec. 3.2-3.3 y Sec. 9.
! Decision metodologica clave: N_e creados ~ Binom(n_F, p_F);
!   extraidos ~ Binom(n_F, p_F*EEE), con EEE=0.99 en Ar.
!=======================================================================
program red100_nest
  use constants
  use mod_tnr_to_e
  use xsections_nest
  use mod_stats          ! poisson_prob, binomial_prob
  use flux
  implicit none

  integer, parameter :: n_T = 500, n_E = 2000, n_ion = 15
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE
  real(dp) :: integrando_kop, integrando_mue, integrando_comb
  integer :: i_T, i_E, i_bin, u_out
  real(dp) :: tasa_Kop, tasa_Mue, tasa_Comb
  real(dp) :: atoms_per_kg, sec_per_day
  real(dp), allocatable :: array_tasa_Comb(:)
  real(dp) :: peso, prob_c, prob_e, p_F
  real(dp) :: tasa_creados(0:n_ion), tasa_extraidos(0:n_ion)
  character(len=250) :: outdir, filename
  integer :: n, n_F
  real(dp) :: dummy
  real(dp) :: R_tot_kop, R_tot_mue, R_tot_comb, sum_cre, sum_ext, mean_cre, mean_ext

  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'

  call inicializar_nest()
  allocate(array_tasa_Comb(n_T))

  ! Carga débil
  QV2 = (-N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge)**2

  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp

  dummy = flujo_diferencial(1.0_dp)

  T_nr_min = 0.20_dp / 1000.0_dp   ! 0.20 keV = suelo de NEST v2.4.0 (antes 0.21)
  T_nr_max = 10.0_dp / 1000.0_dp
  dT = (T_nr_max - T_nr_min) / (n_T - 1)
  dT_keV = dT * 1000.0_dp
  dE = E_nu_max / (n_E - 1)

  ! ================== FASE 1: Espectro continuo ==================
  filename = trim(outdir) // 'espectro_continuoAr.dat'
  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# T_nr[keV]   Tasa_Kop   Tasa_Mue   Tasa_Comb'

  R_tot_kop = 0.0_dp; R_tot_mue = 0.0_dp; R_tot_comb = 0.0_dp

  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Kop = 0.0_dp; tasa_Mue = 0.0_dp; tasa_Comb = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle
        peso = merge(0.5_dp, 1.0_dp, i_E == 1 .or. i_E == n_E)

        integrando_kop = (phi_total * kopeikin_spectrum(E_nu) / spectrum_integral) * dsigma_dT(E_nu, T_nr)
        integrando_mue = (phi_total * (fission_frac(1)*mueller_spectrum(E_nu,1) + &
                                       fission_frac(2)*mueller_spectrum(E_nu,2) + &
                                       fission_frac(3)*mueller_spectrum(E_nu,3) + &
                                       fission_frac(4)*mueller_spectrum(E_nu,4)) &
                          / spectrum_integral) * dsigma_dT(E_nu, T_nr)
        integrando_comb = flujo_diferencial(E_nu) * dsigma_dT(E_nu, T_nr)

        tasa_Kop  = tasa_Kop  + integrando_kop * dE * peso
        tasa_Mue  = tasa_Mue  + integrando_mue * dE * peso
        tasa_Comb = tasa_Comb + integrando_comb * dE * peso
     end do
     tasa_Kop  = tasa_Kop  * atoms_per_kg * sec_per_day
     tasa_Mue  = tasa_Mue  * atoms_per_kg * sec_per_day
     tasa_Comb = tasa_Comb * atoms_per_kg * sec_per_day
     array_tasa_Comb(i_T) = tasa_Comb
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)   ! trapecio en T
     R_tot_kop  = R_tot_kop  + tasa_Kop  * dT_keV * peso
     R_tot_mue  = R_tot_mue  + tasa_Mue  * dT_keV * peso
     R_tot_comb = R_tot_comb + tasa_Comb * dT_keV * peso
     write(u_out, '(F10.4, 3ES15.6)') T_nr*1000.0_dp, tasa_Kop, tasa_Mue, tasa_Comb
  end do
  close(u_out)

  ! ================== FASE 2: Ionización (Poisson + Binomial) ==================
  filename = trim(outdir) // 'ionization_electrones_Ar.dat'
  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# N_e   creados [Evts/kg/dia]   extraidos [Evts/kg/dia]'

  tasa_creados   = 0.0_dp
  tasa_extraidos = 0.0_dp

  ! Fluctuacion de N_e con el modelo de NEST (sub-Poissoniano):
  !   N_e creados   ~ Binomial(n_F, p_F)          p_F = 1 - F(T)
  !   N_e extraidos ~ Binomial(n_F, p_F*EEE)      (adelgazado binomial con EEE)
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     call obtener_nest_binomial(T_nr * 1000.0_dp, n_F, p_F)
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)

     do n = 0, n_ion
        prob_c = binomial_prob(n, n_F, p_F)
        if (prob_c > 0.0_dp) &
           tasa_creados(n) = tasa_creados(n) + array_tasa_Comb(i_T) * dT_keV * peso * prob_c

        prob_e = binomial_prob(n, n_F, p_F * EEE)
        if (prob_e > 0.0_dp) &
           tasa_extraidos(n) = tasa_extraidos(n) + array_tasa_Comb(i_T) * dT_keV * peso * prob_e
     end do
  end do

  do i_bin = 0, n_ion
     write(u_out, '(I5, 2ES15.6)') i_bin, tasa_creados(i_bin), tasa_extraidos(i_bin)
  end do
  close(u_out)

  sum_cre = sum(tasa_creados)
  sum_ext = sum(tasa_extraidos)
  ! sumas pesadas por N_e (numero medio de electrones, no de eventos)
  mean_cre = 0.0_dp; mean_ext = 0.0_dp
  do n = 0, n_ion
     mean_cre = mean_cre + real(n,dp) * tasa_creados(n)
     mean_ext = mean_ext + real(n,dp) * tasa_extraidos(n)
  end do

  write(*,'(/,A)') '=================== VALIDACION red100_nest (vs arXiv:2411.18641) ==================='
  write(*,'(A,ES13.5,A)') ' Tasa CEvNS integrada 0.21-10 keV (Comb) = ', R_tot_comb, ' eventos/(kg dia)'
  write(*,'(A,ES13.5)')   '   idem solo Kopeikin (KI)                = ', R_tot_kop
  write(*,'(A,ES13.5)')   '   idem solo Mueller                     = ', R_tot_mue
  write(*,'(A,F10.4,A,F7.3,A)') ' dR/dT max (espectro continuo)           = ', maxval(array_tasa_Comb), &
       ' eventos/(kg dia keV) en T=', T_nr_min*1000.0_dp + (maxloc(array_tasa_Comb,1)-1)*dT_keV, ' keV'
  write(*,'(A,F8.4)')     '   Comb/Kop (razon de tasas totales)     = ', R_tot_comb/max(R_tot_kop,1.0e-30_dp)
  write(*,'(/,A)')        ' Cierre de la estadistica de N_e (n_ion = 15):'
  write(*,'(A,ES13.5)')   '   Sum tasa_creados(0:15) [eventos]   = ', sum_cre
  write(*,'(A,ES13.5)')   '   Sum tasa_extraidos(0:15) [eventos] = ', sum_ext
  write(*,'(A,F8.4,A)')   '   Sum_cre / R_tot_comb (conservacion de eventos) = ', &
       sum_cre/max(R_tot_comb,1.0e-30_dp), '   (~1 OK; <1 => cola perdida por n_ion)'
  write(*,'(A,F10.4)')    '   <Ne> creados   = Sum n*tasa_creados   = ', mean_cre
  write(*,'(A,F10.4)')    '   <Ne> extraidos = Sum n*tasa_extraidos = ', mean_ext
  write(*,'(A,F8.4,A)')   '   <Ne>_ext / <Ne>_cre = ', mean_ext/max(mean_cre,1.0e-30_dp), &
       '   (deberia ~ EEE = 0.328)'
  write(*,'(/,A)')        ' creados / extraidos por bin (Evts/kg/dia):'
  write(*,'(A)')          '   Ne    creados        extraidos'
  do i_bin = 0, 10
     write(*,'(I5,2ES15.6)') i_bin, tasa_creados(i_bin), tasa_extraidos(i_bin)
  end do
  write(*,'(A,/)')        '==================================================================================='

end program red100_nest
