!=======================================================================
! Parametros fisicos y de RED-100 con xenon liquido (arXiv:2411.18641).
! Los nombres *_Ge vienen del codigo de CONUS+; aqui son del xenon.
!=======================================================================
module constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! ==================== FÍSICA FUNDAMENTAL ====================
  real(dp), parameter :: GF     = 1.1663787d-11    ! Fermi (MeV⁻²)
  real(dp), parameter :: pi     = 3.141592653589793_dp
  real(dp), parameter :: NA     = 6.02214076d23    ! Avogadro (mol⁻¹)
  real(dp), parameter :: hbarc2 = 3.89379d-22      ! (ℏc)² en cm²·MeV²
  real(dp), parameter :: s2w    = 0.23857_dp
  real(dp), parameter :: conv   = 1.0_dp

  ! ==================== BLANCO: XENÓN LÍQUIDO ====================
  real(dp), parameter :: amu  = 931.49410242_dp   ! MeV / u
  real(dp), parameter :: A_Ge = 131.293_dp
  real(dp), parameter :: Z_Ge = 54.0_dp
  ! N promedio de los isotopos naturales (abundancias del NIST)
  real(dp), parameter :: N_Ge = 77.3879_dp
  real(dp), parameter :: M_Ge = A_Ge * amu        ! masa nuclear (MeV)

  ! ==================== DETECTOR RED-100 ====================
  real(dp), parameter :: total_mass_kg   = 126.0_dp   ! masa activa
  real(dp), parameter :: EEE             = 0.328_dp   ! eficiencia de extraccion (32.8 +/- 2.8 %)

  ! Exposicion con el reactor encendido en kg*dia (volumen fiducial)
  real(dp), parameter :: exposure_ON_kgd = 192.0_dp

  ! Fraccion de senal que queda tras los cortes, por bin de N_e reconstruido
  ! (digitalizada de RED-100; en total 25,5 %)
  real(dp), parameter :: eff_ROI(4:7) = (/ 0.1369_dp, 0.3271_dp, 0.6105_dp, 0.7373_dp /)

  ! ==================== FLUJO ====================
  real(dp), parameter :: phi_total     = 1.4d13     ! [nu / cm^2 s] flujo total
  real(dp), parameter :: E_nu_min_flux = 2.0_dp
  real(dp), parameter :: E_nu_max      = 10.0_dp

  real(dp) :: QV2 = 1.0_dp
end module constants
