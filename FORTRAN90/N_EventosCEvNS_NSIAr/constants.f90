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

  ! ------------------------------------------------------------
  ! Parámetros del material: ARGÓN LÍQUIDO (LAr)
  ! *Nota: Se conservan los nombres de variables (ej. A_Ge, Z_Ge) 
  ! para mantener compatibilidad modular con xsections y quenching.
  ! ------------------------------------------------------------
  real(dp), parameter :: A_Ge = 39.948_dp         ! Masa atómica (g/mol) para Ar
  real(dp), parameter :: M_Ge = A_Ge * amu        ! Masa nuclear (MeV)
  integer, parameter :: Z_Ge = 18                 ! Número atómico para Ar
  real(dp), parameter :: N_Ge = 22.0_dp           ! Número de neutrones para Ar

  ! Factor de quenching (modelo de Lindhard)
  real(dp), parameter :: k_lind = 0.1971_dp       ! Parámetro k para Argón
  real(dp), parameter :: alpha = 11.5_dp * Z_Ge**(-7.0_dp/3.0_dp)

  ! Parámetros de resolución energética e ionización
  real(dp), parameter :: sigma0 = 100.0d-6        ! 100 eV
  real(dp), parameter :: eta = 23.6d-6            ! 23.6 eV (energía media por par e- en LAr)
  real(dp), parameter :: Fano = 0.10_dp           ! Factor de Fano para Argón

  ! ------------------------------------------------------------
  ! Flujo de antineutrinos (RED-100: reactor Kalinin a 19 m)
  ! ------------------------------------------------------------
  real(dp), parameter :: phi_total = 1.4d13       ! Flujo total (cm^-2 s^-1)

  ! ------------------------------------------------------------
  ! Exposición y masa del detector (Proyección RED-100 en Argón)
  ! ------------------------------------------------------------
  real(dp), parameter :: exposure_kg_d = 331.0_dp           ! Se mantiene para comparar (kg·días)
  real(dp), parameter :: exposure_kg_s = exposure_kg_d * 86400.0_dp
  real(dp), parameter :: total_mass_kg = 62.0_dp            ! Masa activa proyectada (kg)
  real(dp), parameter :: exposure_atoms_s = (NA / A_Ge) * 1000.0_dp * exposure_kg_s

  ! ------------------------------------------------------------
  ! Binning de energía reconstruida y Umbral
  ! Umbral definido por 4 electrones de ionización
  ! ------------------------------------------------------------
  integer, parameter :: nbins = 19 !conus+
  real(dp), parameter :: Erec_width = 100d-6                
  real(dp), parameter :: EEE = 1.0_dp                     ! Eficiencia Extracción LAr (~100%)
  real(dp), parameter :: threshold = 4.0_dp * eta / EEE   ! Sube automáticamente a ~94.4 eVee
  real(dp), parameter :: Erec_min = threshold

  ! ------------------------------------------------------------
  ! Mallas de integración
  ! ------------------------------------------------------------
  integer, parameter :: npts_Eer = 2000
  real(dp), parameter :: Eer_min = 2.96d-6
  real(dp), parameter :: Eer_max = 0.002_dp                 ! 2 keVee
  integer, parameter :: npts_nu = 2000
  real(dp), parameter :: E_nu_max = 10.0_dp
  real(dp), parameter :: E_nu_min_flux = 2.0_dp
  integer, parameter :: npts_gauss = 1000

  ! Variable global para escalar la carga débil
  real(dp) :: QV2 = 1.0_dp
end module constants
