program int_resolu_all
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846d0

  ! Constantes físicas
  real(dp), parameter :: eta = 2.96d-6          ! MeV (2.96 eV)
  real(dp), parameter :: F_Fano = 0.1096d0
  real(dp), parameter :: sigma0 = 2.038d-5      ! MeV (20.38 eV)
  real(dp), parameter :: E_nu = 10.0d0          ! MeV
  real(dp), parameter :: M = 67663.72d0         ! MeV
  integer, parameter :: A = 76
  integer, parameter :: Z = 32
  real(dp), parameter :: k_lind = 0.162d0

  ! Bins de energía reconstruida
  integer, parameter :: nbins = 19
  real(dp), parameter :: bin_start = 160.0d-6, bin_width = 10.0d-6

  ! Muestreo de TN
  integer, parameter :: npoints = 1000
  real(dp) :: Tmax, dTN, TN, Tmin
  integer :: i, j, k

  real(dp) :: Eer, sres, Gval, sumG, dE
  real(dp) :: Ereco_low, Ereco_high, Kval(nbins)
  integer, parameter :: nE = 1000
  integer :: u_lind, u_noqf
  character(len=256) :: header

  ! Nombres de archivo fijos (los que usa el script Python)
  character(len=*), parameter :: path_lind = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/K_allQF.dat"
  character(len=*), parameter :: path_noqf = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/K_all.dat"

  ! Cálculo de Tmax (energía máxima de retroceso para neutrino de 10 MeV)
  Tmax = 2.0d0 * E_nu**2 / (M + 2.0d0 * E_nu)
  Tmin = 1.0d-4  ! 100 eV
  dTN = (Tmax - Tmin) / real(npoints - 1, dp)

  print *, "Tmax [MeV] =", Tmax, " | sigma0 =", sigma0, " MeV"
  print *, "Generando archivos:"
  print *, "  - ", trim(path_lind), " (con Lindhard)"
  print *, "  - ", trim(path_noqf), " (sin quenching, QF=1)"

  ! Encabezado: primera columna, luego las 19 columnas de los bins
  header = "# Eer[MeV]"
  do j = 1, nbins
     write(header(len_trim(header)+1:), '(A,I2.2)') "   K_bin", j
  end do

  ! Abrir ambos archivos
  open(newunit=u_lind, file=path_lind, status="replace")
  open(newunit=u_noqf, file=path_noqf, status="replace")
  write(u_lind, '(A)') trim(header)
  write(u_noqf, '(A)') trim(header)

  ! Bucle sobre TN
  do i = 0, npoints - 1
     TN = Tmin + i * dTN

     ! --- CASO 1: Con quenching de Lindhard (archivo K_allQF.dat) ---
     Eer = QF_TN(TN) * TN
     sres = sqrt(sigma0**2 + F_Fano * eta * Eer)
     call compute_K(Eer, sres, Kval)
     write(u_lind, '(ES15.8E3,1X,19(ES15.8E3,1X))') Eer, Kval

     ! --- CASO 2: Sin quenching, QF = 1 (archivo K_all.dat) ---
     Eer = TN   ! porque QF=1 => Eer = TN
     sres = sqrt(sigma0**2 + F_Fano * eta * Eer)
     call compute_K(Eer, sres, Kval)
     write(u_noqf, '(ES15.8E3,1X,19(ES15.8E3,1X))') Eer, Kval
  end do

  close(u_lind)
  close(u_noqf)
  print *, "Archivos generados correctamente."

contains

  ! Función de quenching de Lindhard
  function QF_TN(TN) result(QF)
    real(dp), intent(in) :: TN
    real(dp) :: QF, eps, g_eps
    eps = 11.5d0 * Z**(-7.0d0/3.0d0) * (TN * 1.0d3)   ! TN en MeV -> keV
    g_eps = 3.0d0 * eps**0.15d0 + 0.7d0 * eps**0.6d0 + eps
    QF = (k_lind * g_eps) / (1.0d0 + k_lind * g_eps)
  end function QF_TN

  ! Subrutina que integra la gaussiana en los 19 bins para un Eer y sres dados
  subroutine compute_K(Eer, sres, K)
    real(dp), intent(in) :: Eer, sres
    real(dp), intent(out) :: K(:)
    integer :: j, kk
    real(dp) :: Ereco_low, Ereco_high, dE, sumG, Gval, Ereco
    integer, parameter :: nE_local = 1000
    do j = 1, size(K)
       Ereco_low  = bin_start + (j-1) * bin_width
       Ereco_high = Ereco_low + bin_width
       dE = (Ereco_high - Ereco_low) / real(nE_local - 1, dp)
       sumG = 0.0d0
       do kk = 0, nE_local - 1
          Ereco = Ereco_low + kk * dE
          Gval = (1.0d0 / (sqrt(2.0d0 * pi) * sres)) * exp(-0.5d0 * ((Ereco - Eer) / sres)**2)
          if (kk == 0 .or. kk == nE_local - 1) then
             sumG = sumG + 0.5d0 * Gval
          else
             sumG = sumG + Gval
          end if
       end do
       K(j) = sumG * dE
    end do
  end subroutine compute_K

end program int_resolu_all
