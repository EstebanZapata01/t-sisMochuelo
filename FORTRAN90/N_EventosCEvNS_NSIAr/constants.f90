! Parametros para el blanco de argon liquido (Ar-40); clon de la carpeta
! de Xe con los parametros de Ar.
! Fuentes de los parametros de Ar:
!   ref.[46] = D. Akimov et al. (RED-100), Physics 5, 492 (2023),
!              DOI 10.3390/physics5020034 (open access).
!   ReD 2025 = arXiv:2510.16404 (Q_y de NR de Ar medido en 2.4-7.6 keV).
! Cambios respecto a la carpeta de Xe (y su justificacion):
!   - A,Z,N,M = Ar-40 (39.948, 18, 22).
!   - masa activa = 62 kg  (ref.[46], proximo montaje; NO 126).
!   - EEE = 0.99  (ref.[46]: umbral de emision ~0.2 kV/cm en Ar vs ~1.8
!     en Xe -> extraccion ~100 %. ReD: 3.8 kV/cm -> "100 %". DS-50 >99.9 %).
!   - eff_ROI(1:7) = 1  (idealizacion: no hay Fig. 6 de argon).
!   - T_nr_max = 3.5 keV en los chi2 (cinematica E_nu=8 MeV -> 3.44 keV +
!     ref.[46]; trunca cola E_nu 8-10 MeV, <1 %).
!   - ROI en N_e: 1..5 (arXiv:2411.18641 SVII "below five ionization
!     electrons"). INCONSISTENCIA: ref.[46] dice "less than four". Se usa 1..5.
!   - phi_total y el espectro no cambian (mismo reactor KNPP, mismo
!     hibrido Kopeikin+Mueller).
! OJO: el fondo de argon NO se simula: se toma solo el nivel que declara la
!      ref.[46] (S/sqrt(B) ~ 4 a 62 kg*dia) en chi2_bkg_nest.f90. Su composicion
!      (p. ej. apilamiento de electron unico) la deja ABIERTA la ref.[46].
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
  ! Clon del pipeline de Xe (../N_EventosCEvNS_NSIXe) adaptado a Ar-40.
  ! Los nombres siguen siendo *_Ge por herencia (proyecto original de Ge);
  ! aquí contienen valores de ARGÓN. Cambian solo A, Z, N, M respecto a Xe.
  real(dp), parameter :: amu  = 931.49410242_dp   ! MeV / u
  real(dp), parameter :: A_Ge = 39.948_dp         ! Ar natural (dominado por 40Ar)
  real(dp), parameter :: Z_Ge = 18.0_dp
  real(dp), parameter :: N_Ge = 22.0_dp
  real(dp), parameter :: M_Ge = A_Ge * amu        ! masa nuclear (~37218 MeV)

  ! ==================== DETECTOR RED-100 con Ar (arXiv:2411.18641 §VII) ====================
  ! §VII plantea el cambio Xe->Ar como pregunta abierta. Mismo campo de deriva
  ! (218 V/cm, "same electric field strength"); Ar da mayor EEE y menor light
  ! yield de electroluminiscencia. rho(LAr) ~ 1.40 g/cm3 (~87 K).
  real(dp), parameter :: total_mass_kg   = 62.0_dp    ! ref.[46]: masa del proximo montaje de Ar
  real(dp), parameter :: mass_Ar_kg      = 62.0_dp    ! alias explicito para chi2_bkg_nest
  real(dp), parameter :: dias_exposicion = 331.0_dp   ! OBSOLETO: son kg*dia (ver exposure_ON_kgd)
  real(dp), parameter :: livetime_frac   = 0.60_dp    ! informativo

  ! Nivel de fondo de la ref.[46]: la fluctuacion del fondo es ~4x menor que la
  ! senal CEvNS esperada a 62 kg*dia (1 dia), S/sqrt(B) ~ 4. Es el UNICO fondo de Ar.
  real(dp), parameter :: SB_ref46 = 4.0_dp

  ! EEE (eficiencia de extracción liquido->gas) de Ar.
  ! ReD arXiv:2510.16404 §2: "the extraction field guarantees 100% extraction
  ! efficiency". DarkSide-50 mide >99.9% a campos de extracción altos. §VII de
  ! arXiv:2411.18641 dice que Ar extrae MEJOR que el 0.328 de Xe al mismo campo.
  ! Se usa 0.99. Salvedad: depende de un campo de extracción alto que RED-100
  ! con Ar no tiene publicado; para un sistematico, bajar a ~0.90.
  real(dp), parameter :: EEE = 0.99_dp

  ! Exposicion de referencia en kg*dia. Se reusa el 192 (FV) de Xe como
  ! linea base para la comparacion; en la sensibilidad se escanea la exposicion.
  real(dp), parameter :: exposure_ON_kgd = 192.0_dp

  ! Retencion de senal por bin de N_e.
  ! IDEALIZACION COMUN Xe/Ar (metodologia acordada): sin cortes de seleccion,
  ! eff_ROI = 1 en todo el ROI. No existe una Fig. 6 de Ar que digitalizar.
  ! ROI real de Ar (§VII arXiv:2411.18641: "below five ionization electrons",
  ! ref. [46]) = N_e = 1..5, usado por chi2red100_nest.f90 (NE_LO..NE_HI).
  ! Se declara 1..7 para que los programas de diagnostico en PE heredados de Xe
  ! (mainred100_nest.f90 indexa 4..7) sigan compilando; da igual porque todo = 1.
  real(dp), parameter :: eff_ROI(1:7) = 1.0_dp

  ! ==================== FLUJO EXPLICITO (idéntico a Xe) ====================
  ! Mismo espectro híbrido Kopeikin(<2 MeV)+Huber-Mueller(>=2 MeV), misma norma.
  real(dp), parameter :: phi_total     = 1.4d13     ! [nu / cm^2 s] flujo TOTAL
  real(dp), parameter :: E_nu_min_flux = 2.0_dp
  real(dp), parameter :: E_nu_max      = 10.0_dp

  real(dp) :: QV2 = 1.0_dp   ! q_eff^2 (NSI); el programa principal debe asignarlo antes de llamar a dsigma_dT
end module constants
