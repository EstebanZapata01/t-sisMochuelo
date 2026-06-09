program test_integral_mueller
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n_iso = 4
  real(dp), parameter :: fission_frac(n_iso) = [0.58_dp, 0.07_dp, 0.30_dp, 0.05_dp]
  real(dp), parameter :: U235(6)  = [3.217_dp, -3.111_dp, 1.395_dp, -0.369_dp, 0.04445_dp, -0.002053_dp]
  real(dp), parameter :: U238(6)  = [0.4833_dp, 0.1927_dp, -0.1283_dp, -0.006762_dp, 0.002233_dp, -0.0001536_dp]
  real(dp), parameter :: Pu239(6) = [6.413_dp, -7.432_dp, 3.535_dp, -0.882_dp, 0.1025_dp, -0.00455_dp]
  real(dp), parameter :: Pu241(6) = [3.251_dp, -3.204_dp, 1.428_dp, -0.3675_dp, 0.04254_dp, -0.001896_dp]
  real(dp) :: E_nu, integral, h, x
  integer :: i, j, n
  real(dp), external :: espectro_total

  n = 10000
  h = (10.0_dp - 2.0_dp) / n
  integral = 0.0_dp
  do i = 0, n
     x = 2.0_dp + i * h
     if (i == 0 .or. i == n) then
        integral = integral + 0.5_dp * espectro_total(x)
     else
        integral = integral + espectro_total(x)
     end if
  end do
  integral = integral * h
  print *, 'Integral del espectro total de Mueller (2-10 MeV) = ', integral
  print *, 'Valores a 5 MeV:'
  print *, '  espectro_total(5.0) = ', espectro_total(5.0_dp)
  print *, '  (phi_total / 6) = ', 1.5d13/6.0_dp
  print *, '  flujo diferencial con tu normalización = ', (1.5d13/6.0_dp) * espectro_total(5.0_dp)
  print *, '  flujo diferencial con normalización por integral real = ', 1.5d13 * espectro_total(5.0_dp) / integral

contains
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
    if (E_nu < 2.0_dp .or. E_nu > 10.0_dp) then
       rho = 0.0_dp
    else
       rho = 0.0_dp
       do i = 1, n_iso
          rho = rho + fission_frac(i) * mueller_spectrum(E_nu, i)
       end do
    end if
  end function espectro_total
end program test_integral_mueller
