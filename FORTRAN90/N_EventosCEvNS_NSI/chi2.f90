!=======================================================================
! Programa: chi2_sin2theta
!   Calcula χ² para el ángulo de mezcla débil sin²θ_W usando los datos
!   de CONUS+. Incluye minimización sobre el parámetro de nuisance α.
!   Precalcula R_unit (eventos por kg para Q_W = 1) y luego barre sin²θ.
!=======================================================================
program chi2_sin2theta
  use constants
  use quenching
  use xsections
  use flux
  use resolution
  implicit none

  integer :: bin, step, u_out
  integer, parameter :: n_steps = 3000
  real(dp) :: R_exp(nbins), sigma_exp(nbins), R_unit(nbins), R_th(nbins)
  real(dp) :: Eer_vals(npts_Eer), K_matrix(npts_Eer, nbins), S_vals(npts_Eer)
  real(dp) :: s2w, qw, chi2, factor, alpha_best, sum1, sum2
  real(dp) :: s2w_min, s2w_max
  real(dp), parameter :: sigma_alpha = 0.169_dp   ! 16.9%
  character(len=200) :: outdir, filename

  s2w_min = 0.01_dp
  s2w_max = 0.50_dp
  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  filename = trim(outdir) // 'chi2_sin2theta.dat'

  ! Datos experimentales de CONUS+ (excesos en eventos por kg)
  R_exp = [13.27, 38.06, 29.09, 21.37, 12.27, 4.92, 4.29, 0.06, 9.15, &
           2.42, 3.92, -2.55, 0.56, -10.52, -4.29, -8.78, -0.56, -0.93, 13.52]
  sigma_exp = [22.99, 15.70, 12.64, 12.21, 11.52, 11.58, 11.27, 10.77, 11.02, &
               10.28, 10.34, 10.77, 10.34, 10.15, 10.09, 9.90, 9.84, 9.78, 9.59]

  ! ================== 1. Precalcular R_unit (para QV2 = 1) ==================
  call compute_resolution_matrix(npts_Eer, Eer_min, Eer_max, Eer_vals, K_matrix)

  QV2 = 1.0_dp
  do bin = 1, npts_Eer
     call compute_S(Eer_vals(bin), S_vals(bin))
  end do

  do bin = 1, nbins
     R_unit(bin) = integral_trapecio(Eer_vals, K_matrix(:,bin), S_vals, npts_Eer)
     R_unit(bin) = R_unit(bin) * exposure_atoms_s   ! eventos totales para 327 kg·d
     R_unit(bin) = R_unit(bin) / total_mass_kg       ! eventos por kg
  end do

  write(*,*) 'R_unit(1) = ', R_unit(1)

  ! ================== 2. Barrido sobre sin²θ con α ==================
  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# sin2theta_W    chi2'

  do step = 0, n_steps
     s2w = s2w_min + step * (s2w_max - s2w_min) / n_steps
     qw = -N_Ge/2_dp + (1.0_dp - 4.0_dp*s2w)/2_dp * Z_Ge
     factor = qw**2
     R_th(:) = R_unit(:) * factor

     ! Calcular α_best analíticamente
     sum1 = 0.0_dp
     sum2 = 0.0_dp
     do bin = 1, nbins
        sum1 = sum1 + (R_exp(bin) * R_th(bin)) / (sigma_exp(bin)**2)
        sum2 = sum2 + (R_th(bin)**2) / (sigma_exp(bin)**2)
     end do
     alpha_best = (sum1 - sum2) / (sum2 + 1.0_dp/sigma_alpha**2)

     ! Calcular χ² usando α_best
     chi2 = 0.0_dp
     do bin = 1, nbins
        chi2 = chi2 + ((R_exp(bin) - (1.0_dp + alpha_best)*R_th(bin))**2) / (sigma_exp(bin)**2)
     end do
     chi2 = chi2 + (alpha_best / sigma_alpha)**2

     write(u_out, '(F10.6, 2X, F12.4)') s2w, chi2
  end do

  close(u_out)
  write(*,*) 'Archivo generado: ', trim(filename)

contains

  ! ========== Subrutinas auxiliares (compute_S, integrando_S, integral_trapecio) ==========
  subroutine compute_S(Eer, S)
    real(dp), intent(in) :: Eer
    real(dp), intent(out) :: S
    real(dp) :: a, b, T, h, x, sum_trap
    integer :: i
    T = T_from_Eer(Eer)
    if (T <= 0.0_dp) then
       S = 0.0_dp; return
    end if
    a = (T + sqrt(T**2 + 2.0_dp * M_Ge * T)) / 2.0_dp
    if (a >= E_nu_max) then
       S = 0.0_dp; return
    end if
    if (a < E_nu_min_flux) a = E_nu_min_flux
    b = E_nu_max
    h = (b - a) / (npts_nu - 1)
    sum_trap = 0.5_dp * (integrando_S(a, Eer) + integrando_S(b, Eer))
    do i = 2, npts_nu-1
       x = a + (i-1) * h
       sum_trap = sum_trap + integrando_S(x, Eer)
    end do
    S = sum_trap * h
  end subroutine compute_S

  function integrando_S(E_nu, Eer) result(val)
    real(dp), intent(in) :: E_nu, Eer
    real(dp) :: val
    val = flujo_diferencial(E_nu) * dsigma_dEer(E_nu, Eer)
  end function integrando_S

  function integral_trapecio(x, y, z, n) result(intg)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: n
    real(dp) :: intg
    integer :: j
    intg = 0.0_dp
    do j = 1, n-1
       intg = intg + 0.5_dp * (y(j)*z(j) + y(j+1)*z(j+1)) * (x(j+1) - x(j))
    end do
  end function integral_trapecio

end program chi2_sin2theta
