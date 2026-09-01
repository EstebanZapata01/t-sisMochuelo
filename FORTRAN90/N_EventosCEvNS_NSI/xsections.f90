!=======================================================================
! Archivo : xsections.f90   (carpeta CONUS+ / Germanio)
! Rol     : d(sigma)/dEer CE$\nu$NS ya cambiada de variable a energia de
!           ionizacion Eer via el quenching de Lindhard (modulo quenching).
! Tesis   : REFERENCIA; metodologia.tex Sec. 6. En RED-100 la conversion
!           retroceso -> ionizacion la hace NEST (Tnr_to_e.f90), no un
!           factor de quenching analitico.
! Decision metodologica clave: la NSI entra solo por QV2 = q_eff^2, igual
!   que en la version de RED-100.
!=======================================================================
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
