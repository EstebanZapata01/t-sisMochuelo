!=======================================================================
! Programa: mainred100_nest (Asimov Dataset Absoluto - Flujo Explícito)
!   + bloque de VALIDACION para contrastar con arXiv:2411.18641
! Rol     : espectro SM absoluto de eventos CE$\nu$NS por bin de N_e en la
!           ROI 4..7, con exposicion real (192 kg*dia FV).
! Pipeline: integra flux x xsections_nest x binomial(N_e) x eff_ROI ->
!           datos/eventos_sm_xe.dat (consumido por chi2.f90 indirectamente
!           via red100PE) y la tabla de validacion vs Fig. 3 del paper.
! Tesis   : metodologia.tex Sec. 3.3 (Ec. 19) y Sec. 4.
! Decision metodologica clave: Asimov dN_i = R_i; se aplica eff_ROI bin a
!   bin (NO un 0.25 plano); ventana N_e = 4..7 (fondo de SE).
!=======================================================================
program mainred100_nest
  use constants
  use mod_tnr_to_e
  use xsections_nest
  use mod_stats
  use flux, only: flujo_diferencial, E_nu_max, spectrum_integral
  implicit none

  integer, parameter :: n_T = 500, n_E = 2000
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, integrando
  real(dp) :: tasa_Comb, QW_SM, peso
  integer :: i_T, i_E, n_bin, u_out, k
  real(dp) :: R_SM(4:7)
  integer  :: n_F
  real(dp) :: p_F
  character(len=250) :: outdir, filename

  real(dp) :: atoms_per_kg, sec_per_day
  real(dp) :: dummy, Tmax8, tot_roi, lam_c, lam_e, F_here
  real(dp) :: t_ref(5), dRdT_ref(5), tn_ref(6)

  t_ref  = (/ 0.21_dp, 0.30_dp, 0.50_dp, 1.00_dp, 2.00_dp /)   ! keV, para dR/dT
  tn_ref = (/ 0.30_dp, 0.50_dp, 1.00_dp, 2.00_dp, 3.00_dp, 5.00_dp /) ! keV, para NEST

  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  filename = trim(outdir) // 'eventos_sm_xe.dat'

  call inicializar_nest()
  dummy = flujo_diferencial(1.0_dp)   ! dispara init del flujo

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2

  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp

  R_SM = 0.0_dp
  T_nr_min = 0.20_dp / 1000.0_dp; T_nr_max = 3.0_dp / 1000.0_dp   ! 0.20 keV = suelo de NEST v2.4.0
  dT = (T_nr_max - T_nr_min) / (n_T - 1); dT_keV = dT * 1000.0_dp
  dE = E_nu_max / (n_E - 1)

  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Comb = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle
        peso = merge(0.5_dp, 1.0_dp, i_E == 1 .or. i_E == n_E)
        integrando = flujo_diferencial(E_nu) * dsigma_dT(E_nu, T_nr)
        tasa_Comb = tasa_Comb + integrando * dE * peso
     end do
     tasa_Comb = tasa_Comb * atoms_per_kg * sec_per_day
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)

     call obtener_nest_binomial(T_nr * 1000.0_dp, n_F, p_F)
     do n_bin = 4, 7
        R_SM(n_bin) = R_SM(n_bin) + (tasa_Comb * dT_keV * peso * exposure_ON_kgd * &
                      eff_ROI(n_bin) * binomial_prob(n_bin, n_F, p_F * EEE))
     end do
  end do

  ! ---------- dR/dT del espectro continuo en energias de referencia ----------
  do k = 1, size(t_ref)
     T_nr = t_ref(k) / 1000.0_dp
     tasa_Comb = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle
        peso = merge(0.5_dp, 1.0_dp, i_E == 1 .or. i_E == n_E)
        tasa_Comb = tasa_Comb + flujo_diferencial(E_nu) * dsigma_dT(E_nu, T_nr) * dE * peso
     end do
     dRdT_ref(k) = tasa_Comb * atoms_per_kg * sec_per_day
  end do

  ! T_max = 2 E_nu^2 / (M + 2 E_nu), con E_nu=8 MeV y M_Ge en MeV -> pasar a keV
  Tmax8   = 1000.0_dp * (2.0_dp * 8.0_dp**2) / (M_Ge + 2.0_dp*8.0_dp)
  tot_roi = R_SM(4) + R_SM(5) + R_SM(6) + R_SM(7)

  ! ====================================================================
  !                    BLOQUE DE VALIDACION
  ! ====================================================================
  write(*,'(/,A)') '=================== VALIDACION (vs arXiv:2411.18641) ==================='
  write(*,'(A,F12.5,A)')  ' Q_w = N - (1-4 s2w) Z            = ', QW_SM, '   (N=77, Z=54, Ec.2 con 1/2)'
  write(*,'(A,F12.3)')    ' Q_w^2 (en dsigma_dT)             = ', QW_SM**2
  write(*,'(A,ES12.4,A)') ' Flujo total phi                 = ', phi_total, ' nu/cm2/s'
  write(*,'(A,F10.4,A)')  ' Integral espectro (~nubar/fis)  = ', spectrum_integral, '   (paper: 6.75)'
  write(*,'(A,ES12.5,A)') ' Atomos de Xe por kg             = ', atoms_per_kg, ' /kg'
  write(*,'(A,F8.4,A)')   ' T_max retroceso (E_nu = 8 MeV)  = ', Tmax8, ' keV   (paper: ~1 keV)'
  write(*,'(A,F8.1,A)')   ' Exposicion ON                   = ', exposure_ON_kgd, ' kg*dia   (331 activo / 192 FV)'
  write(*,'(A,4F7.3)')    ' eff_ROI(4:7)                    = ', eff_ROI(4), eff_ROI(5), eff_ROI(6), eff_ROI(7)

  write(*,'(/,A)') ' dR/dT espectro continuo [eventos/(kg dia keV)]  (comparar Fig. 2, curva KI/SM2018):'
  do k = 1, size(t_ref)
     write(*,'(A,F5.2,A,ES13.5)') '   T = ', t_ref(k), ' keV : ', dRdT_ref(k)
  end do

  write(*,'(/,A)') ' NEST: rendimiento y fluctuacion  (comparar con GetQuanta de nestpy):'
  write(*,'(A)')   '   T[keV]   <Ne>creados   F=var/media   p_F     n_F   <Ne>extraidos'
  do k = 1, size(tn_ref)
     lam_c  = obtener_electrones_creados(tn_ref(k))
     F_here = obtener_fano(tn_ref(k))
     call obtener_nest_binomial(tn_ref(k), n_F, p_F)
     lam_e  = lam_c * EEE
     write(*,'(F8.2,F13.4,F13.4,F9.4,I7,F13.4)') tn_ref(k), lam_c, F_here, p_F, n_F, lam_e
  end do

  write(*,'(/,A)') ' Asimov SM en la ROI (con eff_ROI, exposicion y Binomial de NEST):'
  do n_bin = 4, 7
     write(*,'(A,I2,A,ES13.5,A,F6.3,A)') '   Bin Ne=', n_bin, ' : ', R_SM(n_bin), &
          '  eventos   (eff_ROI=', eff_ROI(n_bin), ')'
  end do
  write(*,'(A,ES13.5)') '   TOTAL ROI (Ne 4-7)            = ', tot_roi
  write(*,'(A,/)') '======================================================================='

  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# Bin(Ne)   Eventos_SM_Asimov_Absolutos'
  do n_bin = 4, 7
     write(u_out, '(I5, E15.6)') n_bin, R_SM(n_bin)
  end do
  close(u_out)
  write(*,'(A)') ' -> datos/eventos_sm_xe.dat escrito.'
end program mainred100_nest
