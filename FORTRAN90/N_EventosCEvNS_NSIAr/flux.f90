!=======================================================================
! Módulo flux: espectro de antineutrinos de reactor combinado.
!   - Para E_nu < 2 MeV: Kopeikin (2012) con interpolación lineal (50 puntos).
!   - Para E_nu >= 2 MeV: Mueller et al. (2011).
!   La normalización se hace integrando el espectro combinado.
!=======================================================================
module flux
  use constants, only: dp, phi_total, E_nu_min_flux, E_nu_max
  implicit none

  ! ==================== DATOS DE KOPEIKIN (Tabla 3, 50 puntos) ====================
  integer, parameter :: n_kop = 50
  real(dp), parameter :: E_kop(n_kop) = (/ &
       0.010d0, 0.020d0, 0.035d0, 0.040d0, 0.070d0, 0.100d0, 0.130d0, 0.160d0, &
       0.165d0, 0.180d0, 0.215d0, 0.230d0, 0.280d0, 0.330d0, 0.335d0, 0.350d0, &
       0.390d0, 0.400d0, 0.435d0, 0.440d0, 0.500d0, 0.700d0, 0.900d0, 1.000d0, &
       1.185d0, 1.190d0, 1.250d0, 1.300d0, 1.500d0, 1.700d0, 1.800d0, 1.900d0, &
       2.000d0, 2.250d0, 2.500d0, 2.750d0, 3.000d0, 3.250d0, 3.500d0, 4.000d0, &
       4.500d0, 5.000d0, 5.500d0, 6.000d0, 6.500d0, 7.000d0, 7.500d0, 8.000d0, &
       8.500d0, 9.000d0 /)
  real(dp), parameter :: rho_kop(n_kop) = (/ &
       0.774d-1, 0.301d0, 0.848d0, 0.354d0, 0.989d0, 0.181d1, 0.274d1, 0.385d1, &
       0.326d1, 0.374d1, 0.466d1, 0.408d1, 0.503d1, 0.597d1, 0.430d1, 0.408d1, &
       0.442d1, 0.404d1, 0.439d1, 0.284d1, 0.303d1, 0.314d1, 0.323d1, 0.311d1, &
       0.280d1, 0.224d1, 0.206d1, 0.177d1, 0.159d1, 0.148d1, 0.142d1, 0.136d1, &
       0.130d1, 0.108d1, 0.882d0, 0.733d0, 0.611d0, 0.505d0, 0.411d0, 0.261d0, &
       0.156d0, 0.928d-1, 0.549d-1, 0.317d-1, 0.174d-1, 0.897d-2, 0.366d-2, 0.119d-2, &
       0.299d-3, 0.831d-4 /)

  ! ==================== DATOS DE MUELLER ====================
  integer, parameter :: n_iso = 4
  real(dp), parameter :: fission_frac(n_iso) = (/ 0.717d0, 0.068d0, 0.184d0, 0.031d0 /)
  real(dp), parameter :: U235(6)  = (/ 3.217d0, -3.111d0, 1.395d0, -0.369d0, 0.04445d0, -0.002053d0 /)
  real(dp), parameter :: U238(6)  = (/ 0.4833d0, 0.1927d0, -0.1283d0, -0.006762d0, 0.002233d0, -0.0001536d0 /)
  real(dp), parameter :: Pu239(6) = (/ 6.413d0, -7.432d0, 3.535d0, -0.882d0, 0.1025d0, -0.00455d0 /)
  real(dp), parameter :: Pu241(6) = (/ 3.251d0, -3.204d0, 1.428d0, -0.3675d0, 0.04254d0, -0.001896d0 /)

  ! Variables para normalización
  real(dp) :: spectrum_integral = 0.0d0
  logical :: norm_initialized = .false.

contains

  function kopeikin_spectrum(E_nu) result(rho)
    real(dp), intent(in) :: E_nu
    real(dp) :: rho
    integer :: i
    if (E_nu < E_kop(1) .or. E_nu > E_kop(n_kop)) then
       rho = 0.0d0
       return
    end if
    do i = 1, n_kop-1
       if (E_nu >= E_kop(i) .and. E_nu <= E_kop(i+1)) then
          rho = rho_kop(i) + (rho_kop(i+1)-rho_kop(i)) * (E_nu - E_kop(i)) / (E_kop(i+1)-E_kop(i))
          return
       end if
    end do
    rho = rho_kop(n_kop)
  end function kopeikin_spectrum

  function mueller_spectrum(E_nu, iso) result(s)
    real(dp), intent(in) :: E_nu
    integer, intent(in) :: iso
    real(dp) :: s, log_s
    select case(iso)
    case(1)
       log_s = U235(1) + U235(2)*E_nu + U235(3)*E_nu**2 + U235(4)*E_nu**3 &
             + U235(5)*E_nu**4 + U235(6)*E_nu**5
    case(2)
       log_s = U238(1) + U238(2)*E_nu + U238(3)*E_nu**2 + U238(4)*E_nu**3 &
             + U238(5)*E_nu**4 + U238(6)*E_nu**5
    case(3)
       log_s = Pu239(1) + Pu239(2)*E_nu + Pu239(3)*E_nu**2 + Pu239(4)*E_nu**3 &
             + Pu239(5)*E_nu**4 + Pu239(6)*E_nu**5
    case(4)
       log_s = Pu241(1) + Pu241(2)*E_nu + Pu241(3)*E_nu**2 + Pu241(4)*E_nu**3 &
             + Pu241(5)*E_nu**4 + Pu241(6)*E_nu**5
    end select
    s = exp(log_s)
  end function mueller_spectrum

  function espectro_total(E_nu) result(rho)
    real(dp), intent(in) :: E_nu
    real(dp) :: rho
    integer :: i
    if (E_nu < 0.0d0 .or. E_nu > E_nu_max) then
       rho = 0.0d0
       return
    end if
    if (E_nu < E_nu_min_flux) then
       rho = kopeikin_spectrum(E_nu)
    else
       rho = 0.0d0
       do i = 1, n_iso
          rho = rho + fission_frac(i) * mueller_spectrum(E_nu, i)
       end do
    end if
  end function espectro_total

  function compute_spectrum_integral() result(intg)
    real(dp) :: intg
    integer, parameter :: n_int = 1000
    real(dp) :: dE, E
    integer :: i
    dE = E_nu_max / (n_int - 1)
    intg = 0.5d0 * (espectro_total(0.0d0) + espectro_total(E_nu_max))
    do i = 2, n_int-1
       E = (i-1) * dE
       intg = intg + espectro_total(E)
    end do
    intg = intg * dE
  end function compute_spectrum_integral

  function flujo_diferencial(E_nu) result(flujo)
    real(dp), intent(in) :: E_nu
    real(dp) :: flujo
    if (.not. norm_initialized) then
       spectrum_integral = compute_spectrum_integral()
       norm_initialized = .true.
       write(*,*) '=== DIAGNÓSTICO DEL ESPECTRO ==='
       write(*,*) 'Modo: combinado (Kopeikin <2 MeV, Mueller >=2 MeV)'
       write(*,*) 'Integral del espectro (0-10 MeV): ', spectrum_integral
       write(*,*) '================================'
    end if
    flujo = phi_total * espectro_total(E_nu) / spectrum_integral
  end function flujo_diferencial

end module flux
