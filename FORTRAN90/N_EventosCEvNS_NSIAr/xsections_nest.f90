!=======================================================================
! Archivo : xsections_nest.f90
! Rol     : seccion eficaz diferencial CE$\nu$NS  d(sigma)/dT.
! Pipeline: etapa "seccion eficaz" -> la integran mainred100_nest,
!           red100_nest, red100PE, chi2_ideal_nest.
! Tesis   : metodologia.tex Sec. 2.1, Ecs. (1) y (5).
! Decision metodologica clave:
!   - forma  (GF^2 M / pi) * QV2 * (1 - M T/2Enu^2 - T/Enu).
!   - se OMITE el termino subdominante +T^2/(2 Enu^2)  (~1e-4 en la ROI).
!   - por defecto F^2(q^2) = 1. Si la variable de entorno USE_HELM esta a
!     "1"/"true", se incluye el factor de forma de Helm
!     F(q^2) = 3 j1(qR)/(qR) * exp(-(q s)^2 / 2), R = 1.2 A^{1/3} fm,
!     s = 0.9 fm, q = sqrt(2 M T). Sirve para medir el impacto de la
!     aproximacion F^2 = 1 (analisis de sensibilidad).
!   - la NSI entra SOLO por QV2 = q_eff^2 (lo fijan los programas chi2*).
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
    real(dp) :: F2, q, R, x, qs, j1x
    q  = sqrt(2.0_dp * M_Ge * T)                     ! [MeV]
    R  = 1.2_dp * A_Ge**(1.0_dp/3.0_dp)              ! [fm]
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
       use_helm = (trim(env) == '1' .or. trim(env) == 'true' .or. trim(env) == 'TRUE')
       helm_init = .true.
    end if

    ! Prefactor con 1/π (tu definición)
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

end module xsections_nest
