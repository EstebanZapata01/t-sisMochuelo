module xsections
  use constants, only: dp, GF, pi, M_Ge, QV2, conv
  use quenching, only: QF, dQF_dT, T_from_Eer
  implicit none
contains

  function dsigma_dT(E_nu, T) result(dsdT)
    real(dp), intent(in) :: E_nu, T
    real(dp) :: dsdT, prefactor
    ! CORRECCIÓN FÍSICA: Factor de 4 en el denominador
    prefactor = (GF**2 * M_Ge * QV2) / ( pi)
    if (T <= 0.0_dp .or. T >= 2.0_dp*E_nu**2/(M_Ge+2.0_dp*E_nu)) then
       dsdT = 0.0_dp
    else
       dsdT = prefactor * (1.0_dp - (M_Ge*T)/(2.0_dp*E_nu**2) - T/E_nu)
    end if
  end function dsigma_dT

  function dsigma_dEer(E_nu, Eer) result(dsdEer)
    real(dp), intent(in) :: E_nu, Eer
    real(dp) :: dsdEer, T, dqdT
    T = T_from_Eer(Eer)
    if (T <= 0.0_dp) then
       dsdEer = 0.0_dp
       return
    end if
    dqdT = dQF_dT(T) * T + QF(T)
    dsdEer = dsigma_dT(E_nu, T) * (1.0_dp / dqdT) * conv
  end function dsigma_dEer
end module xsections
