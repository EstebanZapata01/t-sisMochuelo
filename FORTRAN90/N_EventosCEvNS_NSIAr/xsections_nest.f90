!=======================================================================
! Seccion eficaz diferencial de CEvNS,
!   d(sigma)/dT = (GF^2 M / pi) QV2 F^2(q) (1 - M T/2Enu^2 - T/Enu).
! Factor de forma de Helm (Lewin y Smith 1996) activo por defecto;
! USE_HELM=0 usa F^2 = 1. La NSI entra solo por QV2 = q_eff^2.
!=======================================================================
module xsections_nest
  use constants, only: dp, GF, pi, M_Ge, QV2, hbarc2, A_Ge
  implicit none

  ! Factor de conversión de MeV⁻³ a cm²/keV
  real(dp), parameter :: conv_cs = hbarc2 / 1000.0_dp
  real(dp), parameter :: hbarc_MeVfm = 197.3269804_dp
  real(dp), parameter :: s_helm_fm   = 0.9_dp

  logical, save :: helm_init = .false.
  logical, save :: use_helm  = .false.

contains

  function helm_F2(T) result(F2)
    real(dp), intent(in) :: T          ! energia de retroceso [MeV]
    real(dp) :: F2, q, c, a, R, x, qs, j1x
    q  = sqrt(2.0_dp * M_Ge * T)                     ! [MeV]
    c  = 1.23_dp * A_Ge**(1.0_dp/3.0_dp) - 0.60_dp   ! [fm]
    a  = 0.52_dp                                     ! [fm]
    R  = sqrt(c**2 + (7.0_dp/3.0_dp)*pi**2*a**2 - 5.0_dp*s_helm_fm**2)  ! [fm]
    x  = q * R / hbarc_MeVfm
    qs = q * s_helm_fm / hbarc_MeVfm
    if (x < 1.0e-8_dp) then
       F2 = 1.0_dp
       return
    end if
    j1x = (sin(x) - x*cos(x)) / x**2                 ! j1(x) = (sin x - x cos x)/x^2
    F2  = (3.0_dp * j1x / x * exp(-0.5_dp*qs**2))**2
  end function helm_F2

  function dsigma_dT(E_nu, T) result(dsdT)
    real(dp), intent(in) :: E_nu, T
    real(dp) :: dsdT, prefactor, ff2
    character(len=8) :: env

    if (.not. helm_init) then
       call get_environment_variable('USE_HELM', env)
       use_helm = .not. (trim(env) == '0' .or. trim(env) == 'false' .or. trim(env) == 'FALSE')
       helm_init = .true.
    end if

    ! Prefactor con 1/π
    prefactor = (GF**2 * M_Ge * QV2) / pi

    if (T <= 0.0_dp .or. T >= 2.0_dp*E_nu**2/(M_Ge+2.0_dp*E_nu)) then
       dsdT = 0.0_dp
    else
       ff2 = 1.0_dp
       if (use_helm) ff2 = helm_F2(T)
       dsdT = prefactor * ff2 * &
              (1.0_dp - (M_Ge*T)/(2.0_dp*E_nu**2) - T/E_nu) * conv_cs
    end if
  end function dsigma_dT

  ! Fija use_helm sin la variable de entorno
  subroutine set_use_helm(flag)
    logical, intent(in) :: flag
    use_helm = flag
    helm_init = .true.
  end subroutine set_use_helm

end module xsections_nest
