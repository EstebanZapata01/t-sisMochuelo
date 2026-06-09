module resolution
  use constants, only: dp, pi, sigma0, eta, Fano, nbins, Erec_min, Erec_width, npts_gauss
  implicit none
contains

  function sigma_res(Eer) result(sig)
    real(dp), intent(in) :: Eer
    real(dp) :: sig
    sig = sqrt(sigma0**2 + Fano * eta * Eer)
  end function sigma_res

  subroutine compute_resolution_matrix(npts, Emin, Emax, Eer_vals, K_matrix)
    integer, intent(in) :: npts
    real(dp), intent(in) :: Emin, Emax
    real(dp), intent(out) :: Eer_vals(npts), K_matrix(npts, nbins)
    real(dp) :: dE, Eer, sig, E_low, E_high, h_g, x, G, sum_trap
    integer :: i, j, k
    dE = (Emax - Emin) / (npts - 1)
    do i = 1, npts
       Eer = Emin + (i-1) * dE
       Eer_vals(i) = Eer
       sig = sigma_res(Eer)
       do j = 1, nbins
          E_low  = Erec_min + (j-1) * Erec_width
          E_high = E_low + Erec_width
          h_g = (E_high - E_low) / (npts_gauss - 1)
          sum_trap = 0.0_dp
          do k = 0, npts_gauss-1
             x = E_low + k * h_g
             G = exp(-0.5_dp * ((x - Eer)/sig)**2) / (sqrt(2.0_dp * pi) * sig)
             if (k == 0 .or. k == npts_gauss-1) then
                sum_trap = sum_trap + 0.5_dp * G
             else
                sum_trap = sum_trap + G
             end if
          end do
          K_matrix(i, j) = sum_trap * h_g
       end do
    end do
  end subroutine compute_resolution_matrix
end module resolution
