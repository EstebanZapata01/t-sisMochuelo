! Parametros fisicos y del detector RED-100 para el blanco de xenon liquido
! (los nombres *_Ge son historicos, de un proyecto previo de Germanio; aqui
! contienen los valores de Xe). Lo usan flux, xsections_nest, Tnr_to_e,
! mod_detector y los programas chi2*.
module constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! ==================== FÍSICA FUNDAMENTAL ====================
  real(dp), parameter :: GF     = 1.1663787d-11    ! Fermi (MeV⁻²)
  real(dp), parameter :: pi     = 3.141592653589793_dp
  real(dp), parameter :: NA     = 6.02214076d23    ! Avogadro (mol⁻¹)
  real(dp), parameter :: hbarc2 = 3.89379d-22      ! (ℏc)² en cm²·MeV²
  real(dp), parameter :: s2w    = 0.23857_dp       ! sin²θ_W (RGE a baja energia)
  real(dp), parameter :: conv   = 1.0_dp

  ! ==================== BLANCO: XENÓN LÍQUIDO ====================
  real(dp), parameter :: amu  = 931.49410242_dp   ! MeV / u
  real(dp), parameter :: A_Ge = 131.293_dp
  real(dp), parameter :: Z_Ge = 54.0_dp
  ! N = <A_masico> - Z, con <A_masico>=131.3879 (promedio de NUMEROS MASICOS
  ! ponderado por abundancia NIST, distinto del peso atomico estandar 131.293 u,
  ! que es un promedio de MASAS y difiere del promedio de A por el exceso de
  ! masa por nucleon, casi constante entre isotopos).
  real(dp), parameter :: N_Ge = 77.3879_dp
  real(dp), parameter :: M_Ge = A_Ge * amu        ! masa nuclear coherente con A (~122299 MeV)

  ! ==================== DETECTOR RED-100 (arXiv:2411.18641) ====================
  real(dp), parameter :: total_mass_kg   = 126.0_dp   ! masa activa (informativo)
  real(dp), parameter :: dias_exposicion = 331.0_dp   ! OBSOLETO: son kg*dia (ver exposure_ON_kgd)
  real(dp), parameter :: livetime_frac   = 0.60_dp    ! informativo: cociente livetime/tiempo-real
  real(dp), parameter :: EEE             = 0.328_dp   ! eficiencia de extraccion (32.8 +/- 2.8 %)

  ! Exposicion reactor ON en kg*dia (livetime 2.63 dia * masa). 192 = volumen
  ! fiducial (el que usa el chi2 del paper); 331 = volumen activo.
  real(dp), parameter :: exposure_ON_kgd = 192.0_dp

  ! Retencion de senal CEvNS tras los cortes de seleccion, por bin de N_e
  ! RECONSTRUIDO (PE/27): cociente senal despues/antes de cortes de arXiv:2411.18641
  ! (digitalizada de la figura de residuos del paper). Se aplica al bin
  ! reconstruido (mod_detector: eps_ROI_pe, prob_migracion), no al N_e
  ! verdadero; 0 fuera de N_e = 4..7 (ROI = corte duro de 110 a 189 PE).
  real(dp), parameter :: eff_ROI(4:7) = (/ 0.1369_dp, 0.3271_dp, 0.6105_dp, 0.7373_dp /)

  ! ==================== FLUJO EXPLICITO ====================
  real(dp), parameter :: phi_total     = 1.4d13     ! [nu / cm^2 s] flujo TOTAL (6.75 nubar/fision)
  real(dp), parameter :: E_nu_min_flux = 2.0_dp
  real(dp), parameter :: E_nu_max      = 10.0_dp

  real(dp) :: QV2 = 1.0_dp   ! q_eff^2 (NSI); el programa principal debe asignarlo antes de llamar a dsigma_dT
end module constants
