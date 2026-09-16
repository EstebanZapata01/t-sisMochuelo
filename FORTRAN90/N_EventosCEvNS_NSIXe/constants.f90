!=======================================================================
! Archivo : constants.f90   (carpeta Xe)
! Rol     : parametros fisicos y del detector RED-100 para el blanco de
!           XENON liquido (los nombres *_Ge son historicos; contienen Xe).
! Pipeline: base de TODO el pipeline -> lo usan flux, xsections_nest,
!           Tnr_to_e, mod_detector y los programas chi2*.
! Tesis   : metodologia.tex Sec. 2 (marco fisico) y Sec. 9 (registro de
!           cambios Ar vs Xe).
! Decision metodologica clave:
!   - sin2th_W = 0.23857 (RGE a baja energia); QW por Ec. (3)-(4).
!   - EEE = 0.328 (32.8 +/- 2.8 %, arXiv:2411.18641); entra en el
!     adelgazado binomial de Ne (Ec. 17 de metodologia.tex).
!   - eff_ROI(4:7) digitalizada de la Fig. 6 (retencion global ~25 %).
!   - phi_total = 1.4e13 es el flujo TOTAL (6.75 nubar/fision); la
!     normalizacion del espectro vive en flux.f90.
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
  ! (nombres *_Ge heredados de un proyecto previo de Germanio; contienen Xe)
  real(dp), parameter :: amu  = 931.49410242_dp   ! MeV / u
  real(dp), parameter :: A_Ge = 131.293_dp
  real(dp), parameter :: Z_Ge = 54.0_dp
  ! N = <A> - Z (promedio isotopico ponderado por abundancia natural, no la
  ! N=77 de un solo isotopo puro): Q_W es lineal en N, asi que el promedio
  ! de Q_W sobre isotopos es exactamente -N_prom/2 + ...*Z, N_prom = <A>-Z.
  ! Antes N_Ge=77.0_dp (inconsistente con A_Ge=131.293: 131.293-54=77.293,
  ! no 77); efecto ~0.4% en Q_W, ~0.8% en la tasa (propor. a Q_W^2).
  real(dp), parameter :: N_Ge = A_Ge - Z_Ge
  real(dp), parameter :: M_Ge = A_Ge * amu        ! masa nuclear coherente con A (~122299 MeV)

  ! ==================== DETECTOR RED-100 (arXiv:2411.18641) ====================
  real(dp), parameter :: total_mass_kg   = 126.0_dp   ! masa activa (informativo)
  real(dp), parameter :: dias_exposicion = 331.0_dp   ! OBSOLETO: son kg*dia, no dias (ver exposure_ON_kgd)
  real(dp), parameter :: livetime_frac   = 0.60_dp    ! informativo: cociente livetime/tiempo-real
  real(dp), parameter :: EEE             = 0.328_dp    ! eficiencia de extraccion (32.8 +/- 2.8 %)

  ! Exposicion reactor ON ya en kg*dia (livetime 2.63 dia * masa).
  ! 192 = volumen fiducial (lo que usa el chi2 del paper); 331 = volumen activo.
  real(dp), parameter :: exposure_ON_kgd = 192.0_dp

  ! Retencion de senal CEvNS tras los cortes de seleccion, por bin de N_e.
  ! Digitalizada de la Fig. 6 de arXiv:2411.18641 (crece con N_e; media ~25%,
  ! consistente con el "75% signal loss in ROI" del texto). 0 fuera de N_e=4..7:
  ! la ROI es un corte duro (1-3 e- se descartan por el fondo de SE espontaneos;
  ! >7 e- cae fuera de la ROI de ese analisis).
  real(dp), parameter :: eff_ROI(4:7) = (/ 0.138_dp, 0.330_dp, 0.599_dp, 0.719_dp /)

  ! ==================== FLUJO EXPLICITO ====================
  real(dp), parameter :: phi_total     = 1.4d13     ! [nu / cm^2 s] flujo TOTAL (6.75 nubar/fision)
  real(dp), parameter :: E_nu_min_flux = 2.0_dp
  real(dp), parameter :: E_nu_max      = 10.0_dp

  real(dp) :: QV2 = 1.0_dp
end module constants
