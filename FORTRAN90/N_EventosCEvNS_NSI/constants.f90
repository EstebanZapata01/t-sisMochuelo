module constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! Constantes fundamentales
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: GF = 1.1663787d-11
  real(dp), parameter :: hbarc = 197.3269804d0
  real(dp), parameter :: conv = (hbarc*1.0d-13)**2
  real(dp), parameter :: NA = 6.02214076d23
  real(dp), parameter :: amu = 931.4940954_dp
  real(dp), parameter :: MeV2keV = 1000.0_dp

  ! Germanio
  real(dp), parameter :: A_Ge = 72.6_dp
  real(dp), parameter :: M_Ge = A_Ge * amu
  integer, parameter :: Z_Ge = 32
  real(dp), parameter :: N_Ge = 41.0_dp

  ! Lindhard
  real(dp), parameter :: k_lind = 0.162_dp
  real(dp), parameter :: alpha = 11.5_dp * Z_Ge**(-7.0_dp/3.0_dp)

  ! Resolución CONUS+
  real(dp), parameter :: sigma0 = 20.38d-6
  real(dp), parameter :: eta = 2.96d-6
  real(dp), parameter :: Fano = 0.1096_dp

  ! Flujo (sin cambios)
  real(dp), parameter :: phi_total = 1.5d13


  ! Exposición: usar directamente 327 kg·d
  real(dp), parameter :: exposure_kg_d = 119.0_dp          ! kg·d
  real(dp), parameter :: exposure_kg_s = exposure_kg_d * 86400.0_dp   ! kg·s
  real(dp), parameter :: total_mass_kg = 1.0_dp           ! masa total de los detectores
  real(dp), parameter :: exposure_atoms_s = (NA / A_Ge) * 1000.0_dp * exposure_kg_s   ! átomos·s

  ! Bins experimentales
  integer, parameter :: nbins = 19
  real(dp), parameter :: Erec_min = 160.0d-6
  real(dp), parameter :: Erec_width = 10.0d-6
  real(dp), parameter :: threshold = 160.0d-6

  ! Mallas de integración
  integer, parameter :: npts_Eer = 2000
  real(dp), parameter :: Eer_min = 2.96d-6
  real(dp), parameter :: Eer_max = 0.001_dp
  integer, parameter :: npts_nu = 2000
  real(dp), parameter :: E_nu_max = 10.0_dp
  real(dp), parameter :: E_nu_min_flux = 2.0_dp
  integer, parameter :: npts_gauss = 1000

  ! Variable global para escalonamiento de la Carga Débil (se ajusta en main)
  real(dp) :: QV2 = 1.0_dp
end module constants
