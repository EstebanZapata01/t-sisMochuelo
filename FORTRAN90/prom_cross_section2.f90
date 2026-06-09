program avg_cross_section
  implicit none
  ! usar una única precisión doble de 64 bits
  integer, parameter :: dp = selected_real_kind(15,307)

  integer, parameter :: nE = 200         ! puntos para Enu
  integer, parameter :: n_iso = 4        ! número de isótopos
  integer, parameter :: nT = 80          ! número de puntos en T

  real(dp) :: Enu_min
  real(dp), parameter :: Enu_max = 8.0_dp
  real(dp), parameter :: T_min = 1.0e-4_dp
  real(dp), parameter :: T_max = 1.0e-2_dp
  real(dp), parameter :: pi = 3.141592653589793_dp

  ! constantes en misma precisión dp
  real(dp), parameter :: GF = 1.1663787e-11_dp   ! MeV^-2 (ya convertido)
  real(dp), parameter :: hbarc = 197.3269804_dp  ! MeV·fm
  real(dp), parameter :: conv = (hbarc**2) * 1.0e-26_dp  ! MeV^2·cm^2 -> convierte MeV^-3 a cm^2/MeV según contexto

  integer, parameter :: Z = 32
  integer, parameter :: N = 40
  real(dp), parameter :: mp = 938.272_dp
  real(dp), parameter :: mn = 939.565_dp
  real(dp), parameter :: M = Z*mp + N*mn

  ! arrays principales en dp
  real(dp) :: E_grid(nE), lambda_tot(nE)
  real(dp) :: T_grid(nT), avg_xs(nT)
  real(dp), dimension(6,n_iso) :: coeffs
  real(dp), dimension(n_iso) :: fission_frac

  integer :: i, j, k, l, iunit
  character(len=200) :: filename
  real(dp) :: dE, poly, val, sigma, integral

  ! --- Coeficientes polinomiales (tal cual los tenías) ---
  coeffs(:,1) = (/ 3.217d0, -3.111d0, 1.395d0, -3.690d-1, 4.445d-2, -2.053d-3 /)  ! U235
  coeffs(:,2) = (/ 4.833d-1, 1.927d-1, -1.283d-1, -6.762d-3, 2.233d-3, -1.536d-4 /) ! U238
  coeffs(:,3) = (/ 6.413d0, -7.432d0, 3.535d0, -8.820d-1, 1.025d-1, -4.550d-3 /)    ! Pu239
  coeffs(:,4) = (/ 3.251d0, -3.204d0, 1.428d0, -3.675d-1, 4.254d-2, -1.896d-3 /)    ! Pu241

  ! --- CORRECCIÓN: fracciones en el mismo orden que coeffs(:,i) ---
  fission_frac = (/ 0.58_dp, 0.07_dp, 0.30_dp, 0.05_dp /)

  ! --- Malla en T ---
  do j=1,nT
     T_grid(j) = T_min + (T_max - T_min)*real(j-1,dp)/real(nT-1,dp)
  end do

  ! --- Loop principal en T: reconstruir E_grid desde Enu_min(T), calcular lambda_tot y la integral ---
  do j=1,nT
     integral = 0.0_dp
     lambda_tot(:) = 0.0_dp

     ! umbral cinemático para este T (con masa del núcleo M)
     Enu_min = 0.5_dp * ( T_grid(j) + sqrt( T_grid(j)**2 + 2.0_dp*M*T_grid(j) ) )

     ! si Enu_min >= Enu_max, integral = 0
     if (Enu_min >= Enu_max) then
        avg_xs(j) = 0.0_dp
        cycle
     end if

     ! reconstruir grilla en [Enu_min, Enu_max]
     dE = (Enu_max - Enu_min)/real(nE-1,dp)
     do i=1,nE
        E_grid(i) = Enu_min + dE*real(i-1,dp)
     end do

     ! calcular lambda_tot en esta grilla (suma de isotopos, ponderada por fracciones)
     do k = 1, nE
        lambda_tot(k) = 0.0_dp
        do i = 1, n_iso
           poly = 0.0_dp
           do l = 1, 6
              poly = poly + coeffs(l,i) * E_grid(k)**(l-1)
           end do
           val = fission_frac(i) * exp(poly)
           lambda_tot(k) = lambda_tot(k) + val
        end do
     end do

     ! trapecio sobre la grilla local
     do k=1,nE-1
        sigma = 0.5_dp * ( diff_xs(E_grid(k),T_grid(j)) * lambda_tot(k) + &
                           diff_xs(E_grid(k+1),T_grid(j)) * lambda_tot(k+1) )
        integral = integral + sigma * (E_grid(k+1)-E_grid(k))
     end do
     avg_xs(j) = integral
  end do

  ! imprimir y guardar resultados
  print*, 'Primeros valores avg_xs: ', avg_xs(1:min(8,nT))

  filename = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/promcs.dat'
  open(newunit=iunit, file=filename, status='replace', action='write', form='formatted')
  do j=1,nT
     write(iunit,'(1PE25.16,1X,1PE25.16)') T_grid(j), avg_xs(j)
  end do
  close(iunit)
  print *, 'Datos guardados en ', trim(filename)

contains

  ! diff_xs ahora en la misma precisión dp (usa conv en dp)
  real(dp) function diff_xs(Enu, T)
    implicit none
    real(dp), intent(in) :: Enu, T
    real(dp) :: Qw, prefactor

    Qw = real(N,dp) - (1.0_dp - 4.0_dp*0.231_dp) * real(Z,dp)
    prefactor = (GF**2 / (4.0_dp*pi)) * Qw**2 * M   ! unidades MeV^-3
    diff_xs = prefactor * (1.0_dp - (M*T)/(2.0_dp*Enu**2)) * conv
  end function diff_xs

end program avg_cross_section
