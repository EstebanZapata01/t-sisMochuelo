! Driver minimo para las pruebas de cierre numerico del pipeline: llama a
! las MISMAS rutinas que los programas de produccion (obtener_nest_binomial,
! binomial_prob, dsigma_dT, prob_migracion, eff_ROI, flujo_diferencial) y
! vuelca los resultados a texto, sin reimplementar nada.
! Salidas: datos/cierre_{binomial,dsigma,factorizacion}_Xe.dat
program cierre_pipeline
  use constants
  use mod_tnr_to_e
  use xsections_nest
  use mod_stats, only: binomial_prob
  use mod_detector, only: prob_migracion
  use flux, only: flujo_diferencial, E_nu_max
  implicit none

  integer, parameter :: KMAX = 50
  character(len=250) :: outdir
  integer :: u, i_T, i_E, k, j, n_F
  real(dp) :: T_keV, p_F, pext, dPk, a1, a2, a3, a4, aroi
  real(dp) :: pmf(0:KMAX), suma
  real(dp) :: E_nu, dE, Tmax, dT, T_MeV, tasa, num, sigma_lin
  integer, parameter :: N_T1 = 2000, N_ENU = 300, N_TMID = 4000
  real(dp) :: QW_SM

  outdir = '../../datos/'
  call inicializar_nest()
  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2

  ! ================================================================
  ! (a),(b) delta_P(T) y aceptancia a_i^ROI(T), T en [0.2,3.0] keV
  ! ================================================================
  open(newunit=u, file=trim(outdir)//'cierre_binomial_Xe.dat', status='replace')
  write(u,'(A)') '# T_nr[keV]  delta_P  a(Ne>=1) a(Ne>=2) a(Ne>=3) a(Ne>=4)  a_ROI_real(4-7)'
  do i_T = 1, N_T1
     T_keV = 0.2_dp + (i_T - 1) * (3.0_dp - 0.2_dp) / (N_T1 - 1)
     call obtener_nest_binomial(T_keV, n_F, p_F)
     pext = p_F * EEE
     suma = 0.0_dp; a1=0; a2=0; a3=0; a4=0; aroi=0
     do k = 0, min(n_F, KMAX)
        pmf(k) = binomial_prob(k, n_F, pext)
        suma = suma + pmf(k)
        if (k>=1) a1 = a1 + pmf(k)
        if (k>=2) a2 = a2 + pmf(k)
        if (k>=3) a3 = a3 + pmf(k)
        if (k>=4) a4 = a4 + pmf(k)
        if (k>=1) then
           do j = 4, 7
              aroi = aroi + eff_ROI(j) * prob_migracion(j, k) * pmf(k)
           end do
        end if
     end do
     dPk = abs(1.0_dp - suma)
     write(u,'(7ES16.7)') T_keV, dPk, a1, a2, a3, a4, aroi
  end do
  close(u)

  ! ================================================================
  ! (c) C_sigma(E_nu) - 1, dsigma_dT REAL (QV2=QW_SM^2 fijo aqui)
  ! ================================================================
  open(newunit=u, file=trim(outdir)//'cierre_dsigma_Xe.dat', status='replace')
  write(u,'(A)') '# E_nu[MeV]  C_sigma-1   (num=integral midpoint de dsigma_dT real; sigma_lin=(M QW^2/2pi)*Tmax)'
  do i_E = 1, N_ENU
     E_nu = 0.5_dp + (i_E - 1) * (10.0_dp - 0.5_dp) / (N_ENU - 1)
     Tmax = 2.0_dp * E_nu**2 / (M_Ge + 2.0_dp * E_nu)
     dT = Tmax / N_TMID
     num = 0.0_dp
     do i_T = 1, N_TMID
        T_MeV = (real(i_T,dp) - 0.5_dp) * dT
        num = num + dsigma_dT(E_nu, T_MeV) * dT
     end do
     sigma_lin = (GF**2 * M_Ge * QW_SM**2 / (2.0_dp*pi)) * Tmax * conv_cs
     write(u,'(2ES16.7)') E_nu, num/sigma_lin - 1.0_dp
  end do
  close(u)

  ! ================================================================
  ! cierre de factorizacion: R_k(eps) real, QV2 conmutado, N_e VERDADERO
  ! y ROI reconstruida real (4-7), para A_target en {4.0, 0.3}
  ! ================================================================
  open(newunit=u, file=trim(outdir)//'cierre_factorizacion_Xe.dat', status='replace')
  write(u,'(A)') '# A_target   k(Ne_verdadero, 1..15 o j_ROI=4..7 con flag)   espacio(0=verdadero,1=ROI_real)   R_k'
  call escribe_Rk(u, 1.0_dp, QW_SM)
  call escribe_Rk(u, 4.0_dp, QW_SM)
  call escribe_Rk(u, 0.3_dp, QW_SM)
  close(u)

contains

  subroutine escribe_Rk(u_out, A_target, QW0)
    integer, intent(in) :: u_out
    real(dp), intent(in) :: A_target, QW0
    real(dp) :: Rk_true(1:15), Rk_roi(4:7)
    real(dp) :: T_nr, T_nr_min, T_nr_max, dTl, dTl_keV, dEl, peso, tasaC
    integer :: iT, iE, kk, jj, nF2
    real(dp) :: pF2
    QV2 = A_target * QW0**2
    Rk_true = 0.0_dp; Rk_roi = 0.0_dp
    T_nr_min = 0.20_dp/1000.0_dp; T_nr_max = 3.0_dp/1000.0_dp
    dTl = (T_nr_max - T_nr_min) / (800 - 1); dTl_keV = dTl*1000.0_dp
    dEl = E_nu_max / (2000 - 1)
    do iT = 1, 800
       T_nr = T_nr_min + (iT-1)*dTl
       tasaC = 0.0_dp
       do iE = 1, 2000
          if ((iE-1)*dEl < sqrt(M_Ge*T_nr/2.0_dp)) cycle
          peso = merge(0.5_dp, 1.0_dp, iE==1 .or. iE==2000)
          tasaC = tasaC + flujo_diferencial((iE-1)*dEl) * dsigma_dT((iE-1)*dEl, T_nr) * dEl * peso
       end do
       tasaC = tasaC * (NA/A_Ge) * 1000.0_dp * 86400.0_dp
       peso = merge(0.5_dp, 1.0_dp, iT==1 .or. iT==800)
       call obtener_nest_binomial(T_nr*1000.0_dp, nF2, pF2)
       do kk = 1, min(nF2, 15)
          Rk_true(kk) = Rk_true(kk) + tasaC*dTl_keV*peso*binomial_prob(kk, nF2, pF2*EEE)
       end do
       do kk = 1, nF2
          do jj = 4, 7
             Rk_roi(jj) = Rk_roi(jj) + tasaC*dTl_keV*peso*eff_ROI(jj)*prob_migracion(jj,kk)* &
                          binomial_prob(kk, nF2, pF2*EEE)
          end do
       end do
    end do
    do kk = 1, 15
       write(u_out,'(F8.3,I5,I3,ES25.16)') A_target, kk, 0, Rk_true(kk)
    end do
    do jj = 4, 7
       write(u_out,'(F8.3,I5,I3,ES25.16)') A_target, jj, 1, Rk_roi(jj)
    end do
  end subroutine escribe_Rk

end program cierre_pipeline
