!=======================================================================
! Parametros de RED-100 con argon liquido. Igual que la carpeta del xenon
! salvo el nucleo, la masa (62 kg), la EEE y la eficiencia de los cortes.
! Fuentes: Akimov et al., Physics 5, 492 (2023) y ReD, arXiv:2510.16404.
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

  ! ==================== BLANCO: ARGÓN LÍQUIDO ====================
  real(dp), parameter :: amu  = 931.49410242_dp   ! MeV / u
  real(dp), parameter :: A_Ge = 39.948_dp         ! Ar natural (casi todo 40Ar)
  real(dp), parameter :: Z_Ge = 18.0_dp
  real(dp), parameter :: N_Ge = 22.0_dp
  real(dp), parameter :: M_Ge = A_Ge * amu        ! masa nuclear (MeV)

  ! ==================== DETECTOR RED-100 con Ar ====================
  real(dp), parameter :: total_mass_kg   = 62.0_dp    ! masa del montaje con argon
  real(dp), parameter :: mass_Ar_kg      = 62.0_dp

  ! Fondo declarado por RED-100: S/sqrt(B) ~ 4 a 62 kg*dia
  real(dp), parameter :: SB_ref46 = 4.0_dp

  ! Eficiencia de extraccion (ReD y DarkSide-50 miden ~100 %)
  real(dp), parameter :: EEE = 0.99_dp

  ! Exposicion de referencia en kg*dia, la misma del xenon
  real(dp), parameter :: exposure_ON_kgd = 192.0_dp

  ! Sin cortes de seleccion: eficiencia 1 en todos los bins
  real(dp), parameter :: eff_ROI(1:7) = 1.0_dp

  ! ==================== FLUJO (igual que en el xenon) ====================
  real(dp), parameter :: phi_total     = 1.4d13     ! [nu / cm^2 s] flujo total
  real(dp), parameter :: E_nu_min_flux = 2.0_dp
  real(dp), parameter :: E_nu_max      = 10.0_dp

  real(dp) :: QV2 = 1.0_dp
end module constants
