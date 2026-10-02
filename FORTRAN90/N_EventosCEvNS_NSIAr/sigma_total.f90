! Seccion eficaz total sigma(E_nu) para Ar; gemelo de
! N_EventosCEvNS_NSIXe/sigma_total.f90.
program sigma_total
  use constants
  use xsections_nest
  implicit none
  integer, parameter :: N_ENU = 200, N_TMID = 4000
  real(dp) :: QW_SM, E_nu, Tmax, dT, T_MeV, sigma, sigma_lin
  integer :: i_E, i_T, u
  logical :: helm
  character(len=8) :: env
  character(len=300) :: outdir, fname

  outdir = '../../datos/'
  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2

  call get_environment_variable('USE_HELM', env)
  helm = (trim(env) == '1' .or. trim(env) == 'true')
  fname = trim(outdir) // merge('sigma_total_helm_Ar.dat', 'sigma_total_Ar.dat     ', helm)

  open(newunit=u, file=trim(fname), status='replace')
  write(u,'(A)') '# E_nu[MeV]  sigma[cm^2]  sigma_lin_approx[cm^2]'
  do i_E = 1, N_ENU
     E_nu = 0.2_dp + (i_E - 1) * (10.0_dp - 0.2_dp) / (N_ENU - 1)
     Tmax = 2.0_dp * E_nu**2 / (M_Ge + 2.0_dp * E_nu)
     dT = Tmax / N_TMID
     sigma = 0.0_dp
     do i_T = 1, N_TMID
        T_MeV = (real(i_T,dp) - 0.5_dp) * dT
        sigma = sigma + dsigma_dT(E_nu, T_MeV) * dT
     end do
     sigma_lin = (GF**2 / pi) * QW_SM**2 * E_nu**2 * conv_cs
     write(u,'(3ES18.9)') E_nu, sigma, sigma_lin
  end do
  close(u)
  write(*,'(A,A)') ' -> ', trim(fname)
  write(*,'(A,L1,A,ES12.4)') ' USE_HELM=', helm, '   QW_SM^2=', QW_SM**2
end program sigma_total
