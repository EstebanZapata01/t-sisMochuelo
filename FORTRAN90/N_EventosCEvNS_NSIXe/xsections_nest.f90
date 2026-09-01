!=======================================================================
! Archivo : xsections_nest.f90
! Rol     : seccion eficaz diferencial CE$\nu$NS  d(sigma)/dT.
! Pipeline: etapa "seccion eficaz" -> la integran mainred100_nest,
!           red100_nest, red100PE, chi2_ideal_nest.
! Tesis   : metodologia.tex Sec. 2.1, Ecs. (1) y (5).
! Decision metodologica clave:
!   - forma  (GF^2 M / pi) * QV2 * (1 - M T/2Enu^2 - T/Enu).
!   - se OMITE el termino subdominante +T^2/(2 Enu^2)  (~1e-4 en la ROI).
!   - se fija F^2(q^2) = 1  (factor de Helm): qR << 1 en la ROI, y el
!     error se cancela en la comparacion relativa Xe vs Ar.
!   - la NSI entra SOLO por QV2 = q_eff^2 (lo fijan los programas chi2*).
!=======================================================================
module xsections_nest
  use constants, only: dp, GF, pi, M_Ge, QV2, hbarc2
  implicit none

  ! Factor de conversión de MeV⁻³ a cm²/keV
  real(dp), parameter :: conv_cs = hbarc2 / 1000.0_dp

contains

  function dsigma_dT(E_nu, T) result(dsdT)
    real(dp), intent(in) :: E_nu, T
    real(dp) :: dsdT, prefactor

    ! Prefactor con 1/π (tu definición)
    prefactor = (GF**2 * M_Ge * QV2) / pi

    if (T <= 0.0_dp .or. T >= 2.0_dp*E_nu**2/(M_Ge+2.0_dp*E_nu)) then
       dsdT = 0.0_dp
    else
       dsdT = prefactor * (1.0_dp - (M_Ge*T)/(2.0_dp*E_nu**2) - T/E_nu) * conv_cs
    end if
  end function dsigma_dT

end module xsections_nest
