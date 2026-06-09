program spectrum_antineutrinos
  implicit none
  integer, parameter :: n_iso = 4, n_coef = 6, n_points = 400
  real(8), parameter :: Emin = 2.0d0, Emax = 8.0d0
  real(8), parameter :: fission_frac(n_iso) = (/ 0.58d0, 0.07d0, 0.30d0, 0.05d0 /)
  real(8), dimension(n_points)::lambda_tot
  real(8) :: coeff(n_iso, n_coef)
  real(8) :: E, dE, val, poly
  integer :: i, j, k, unit
  character(len=200) :: filename

 
  character(len=5), dimension(n_iso) :: isotopes
  isotopes = (/ 'U235 ', 'U238 ', 'Pu239', 'Pu241' /)

  coeff = reshape( (/ &
     ! U235
     3.217d0, -3.111d0, 1.395d0, -3.690d-1, 4.445d-2, -2.053d-3, &
     ! U238
     4.833d-1, 1.927d-1, -1.283d-1, -6.762d-3, 2.233d-3, -1.536d-4, &
     ! Pu239
     6.413d0, -7.432d0, 3.535d0, -8.820d-1, 1.025d-1, -4.550d-3, &
     ! Pu241
     3.251d0, -3.204d0, 1.428d0, -3.675d-1, 4.254d-2, -1.896d-3  &
     /), (/ n_iso, n_coef /), order=(/2,1/) )

  dE = (Emax - Emin) / real(n_points-1,8)

  do i = 1, n_iso
     filename = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/spectrum_"//trim(isotopes(i))//".dat"
     open(newunit=unit, file=filename, status="replace", action="write")
     write(unit,'(A)') "# E_nu(MeV)   lambda(E) [MeV^-1 per fission]"

     do j = 0, n_points-1
        E = Emin + dE*real(j,8)
        poly = 0.0d0
        do k = 1, n_coef
           poly = poly + coeff(i,k) * E**(k-1)
        end do
        val = fission_frac(i)*exp(poly)
        lambda_tot(j+1)=lambda_tot(j+1)+val
        print*,i
        print*,lambda_tot(j+1),val
        write(unit,'(F8.3,1X,E18.10)') E, val
     end do

     close(unit)
     print *, "Archivo generado:", trim(filename)
  end do
  print*,lambda_tot
end program spectrum_antineutrinos
