module quenching
  use constants, only: dp, k_lind, alpha, MeV2keV
  implicit none
contains

  function QF(T) result(q)
    real(dp), intent(in) :: T
    real(dp) :: q, T_keV, eps, g
    T_keV = T * MeV2keV
    eps = alpha * T_keV
    ! Aproximación de Robinson para g(epsilon)
    g = 3.0_dp * eps**0.15_dp + 0.7_dp * eps**0.6_dp + eps
    q = (k_lind * g) / (1.0_dp + k_lind * g)
  end function QF

  function dQF_dT(T) result(dq)
    real(dp), intent(in) :: T
    real(dp) :: dq, T_keV, eps, g, dg
    T_keV = T * MeV2keV
    eps = alpha * T_keV
    g = 3.0_dp * eps**0.15_dp + 0.7_dp * eps**0.6_dp + eps
    dg = 0.45_dp * eps**(-0.85_dp) + 0.42_dp * eps**(-0.4_dp) + 1.0_dp
    dq = (k_lind * dg * alpha / (1.0_dp + k_lind * g)**2) * MeV2keV
  end function dQF_dT

  function T_from_Eer(Eer) result(T)
    real(dp), intent(in) :: Eer
    real(dp) :: T, T_low, T_high, E_test
    integer :: iter
    
    ! Búsqueda por bisección pura (sin impresiones lentas)
    T_low = 0.0_dp
    T_high = Eer * 20.0_dp 
    
    do iter = 1, 100
       T = (T_low + T_high) / 2.0_dp
       E_test = T * QF(T)
       if (abs(E_test - Eer) < 1.0d-12) exit
       if (E_test < Eer) then
          T_low = T
       else
          T_high = T
       end if
    end do
  end function T_from_Eer

end module quenching
