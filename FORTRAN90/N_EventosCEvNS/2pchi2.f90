!=======================================================================
! Programa: chi2_nsi_2D (Versión CONUS+ / Germanio)
!   - Usa matriz de resolución gaussiana para Germanio[cite: 6].
!   - Datos experimentales de excesos de CONUS+[cite: 2].
!   - Eje X (eps_x): Primera columna, usualmente sabor electrónico (ee).
!   - Eje Y (eps_y): Segunda columna, sabores mu o tau.
!=======================================================================
program chi2_nsi_2D
  use constants
  use quenching
  use xsections
  use flux
  use resolution
  implicit none

  integer :: i, j, bin, u_out, u_conf
  integer, parameter :: n_y = 1000, n_x = 1000
  real(dp) :: R_exp(nbins), sigma_exp(nbins), R_unit(nbins), R_th(nbins)
  real(dp) :: Eer_vals(npts_Eer), K_matrix(npts_Eer, nbins), S_vals(npts_Eer)
  real(dp) :: eps_y, eps_x, chi2, alpha_best, sum1, sum2
  real(dp) :: eps_y_min, eps_y_max, eps_x_min, eps_x_max, deps_y, deps_x
  real(dp), parameter :: sin2th_SM = 0.23857_dp
  real(dp), parameter :: sigma_alpha = 0.169_dp
  real(dp) :: QW_SM, q_nsi_ee, q_nsi_emu, q_nsi_etau, q_eff2
  
  integer :: ipar ! Selector de caso (1-15)
  character(len=200) :: outdir, filename, configfile
  character(len=50)  :: xlabel, ylabel

  ! ================== CONFIGURACIÓN DEL BARRIDO ==================
  ipar = 5

  eps_y_min = -1.0_dp; eps_y_max = 1.0_dp
  eps_x_min = -1.0_dp; eps_x_max = 1.0_dp
  deps_y = (eps_y_max - eps_y_min) / (n_y - 1)
  deps_x = (eps_x_max - eps_x_min) / (n_x - 1)

  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  filename = trim(outdir) // 'chi2_nsi_2Dconus.dat'
  configfile = trim(outdir) // 'nsi_config_conus.txt'

  ! Datos experimentales CONUS+ (Excesos)[cite: 2]
  R_exp = [13.27, 38.06, 29.09, 21.37, 12.27, 4.92, 4.29, 0.06, 9.15, &
           2.42, 3.92, -2.55, 0.56, -10.52, -4.29, -8.78, -0.56, -0.93, 13.52]
  sigma_exp = [22.99, 15.70, 12.64, 12.21, 11.52, 11.58, 11.27, 10.77, 11.02, &
               10.28, 10.34, 10.77, 10.34, 10.15, 10.09, 9.90, 9.84, 9.78, 9.59]

  ! ================== 1. PRECALCULAR R_unit (QV2 = 1) ==================
  ! Para Germanio es obligatorio usar la matriz de resolución[cite: 6]
  call compute_resolution_matrix(npts_Eer, Eer_min, Eer_max, Eer_vals, K_matrix)
  QV2 = 1.0_dp 
  do i = 1, npts_Eer
     call compute_S(Eer_vals(i), S_vals(i))
  end do
  do bin = 1, nbins
     R_unit(bin) = integral_trapecio(Eer_vals, K_matrix(:,bin), S_vals, npts_Eer)
     R_unit(bin) = R_unit(bin) * exposure_atoms_s / total_mass_kg
  end do

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*sin2th_SM)/2.0_dp * Z_Ge

  ! ================== 2. DEFINICIÓN DE CASOS (EE EN EJE X) ==================
  select case(ipar)
  case(1);  xlabel='$\epsilon_{ee}^{dV}$'; ylabel='$\epsilon_{ee}^{uV}$'
  case(2);  xlabel='$\epsilon_{e\mu}^{dV}$'; ylabel='$\epsilon_{e\mu}^{uV}$'
  case(3);  xlabel='$\epsilon_{e\tau}^{dV}$'; ylabel='$\epsilon_{e\tau}^{uV}$'
  case(4);  xlabel='$\epsilon_{ee}^{uV}$'; ylabel='$\epsilon_{e\mu}^{uV}$'
  case(5);  xlabel='$\epsilon_{ee}^{dV}$'; ylabel='$\epsilon_{e\mu}^{dV}$'
  case(6);  xlabel='$\epsilon_{ee}^{uV}$'; ylabel='$\epsilon_{e\mu}^{dV}$'
  case(7);  xlabel='$\epsilon_{ee}^{dV}$'; ylabel='$\epsilon_{e\mu}^{uV}$'
  case(8);  xlabel='$\epsilon_{ee}^{uV}$'; ylabel='$\epsilon_{e\tau}^{uV}$'
  case(9);  xlabel='$\epsilon_{ee}^{dV}$'; ylabel='$\epsilon_{e\tau}^{dV}$'
  case(10); xlabel='$\epsilon_{ee}^{uV}$'; ylabel='$\epsilon_{e\tau}^{dV}$'
  case(11); xlabel='$\epsilon_{ee}^{dV}$'; ylabel='$\epsilon_{e\tau}^{uV}$'
  case(12); xlabel='$\epsilon_{e\mu}^{uV}$'; ylabel='$\epsilon_{e\tau}^{uV}$'
  case(13); xlabel='$\epsilon_{e\mu}^{dV}$'; ylabel='$\epsilon_{e\tau}^{dV}$'
  case(14); xlabel='$\epsilon_{e\mu}^{uV}$'; ylabel='$\epsilon_{e\tau}^{dV}$'
  case(15); xlabel='$\epsilon_{e\mu}^{dV}$'; ylabel='$\epsilon_{e\tau}^{uV}$'
  end select

  open(newunit=u_conf, file=configfile, status='replace')
  write(u_conf, '(A)') trim(xlabel); write(u_conf, '(A)') trim(ylabel); close(u_conf)

  ! ================== 3. BARRIDO 2D CON MINIMIZACIÓN DE ALFA ==================
  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# eps_x   eps_y   chi2'

  do i = 0, n_y-1
     eps_y = eps_y_min + i * deps_y ! Variable Eje Y
     do j = 0, n_x-1
        eps_x = eps_x_min + j * deps_x ! Variable Eje X
        q_nsi_ee = 0.0_dp; q_nsi_emu = 0.0_dp; q_nsi_etau = 0.0_dp

        select case(ipar)
        case(1); q_nsi_ee = (2.0_dp*eps_y + eps_x)*Z_Ge + (eps_y + 2.0_dp*eps_x)*N_Ge 
        case(2); q_nsi_emu = (2.0_dp*eps_y + eps_x)*Z_Ge + (eps_y + 2.0_dp*eps_x)*N_Ge
        case(3); q_nsi_etau = (2.0_dp*eps_y + eps_x)*Z_Ge + (eps_y + 2.0_dp*eps_x)*N_Ge
        case(4); q_nsi_ee = (2.0_dp*eps_x)*Z_Ge + (eps_x)*N_Ge; q_nsi_emu = (2.0_dp*eps_y)*Z_Ge + (eps_y)*N_Ge
        case(5); q_nsi_ee = (eps_x)*Z_Ge + (2.0_dp*eps_x)*N_Ge; q_nsi_emu = (eps_y)*Z_Ge + (2.0_dp*eps_y)*N_Ge
        case(6); q_nsi_ee = (2.0_dp*eps_x)*Z_Ge + (eps_x)*N_Ge; q_nsi_emu = (eps_y)*Z_Ge + (2.0_dp*eps_y)*N_Ge
        case(7); q_nsi_ee = (eps_x)*Z_Ge + (2.0_dp*eps_x)*N_Ge; q_nsi_emu = (2.0_dp*eps_y)*Z_Ge + (eps_y)*N_Ge
        case(8); q_nsi_ee = (2.0_dp*eps_x)*Z_Ge + (eps_x)*N_Ge; q_nsi_etau = (2.0_dp*eps_y)*Z_Ge + (eps_y)*N_Ge
        case(9); q_nsi_ee = (eps_x)*Z_Ge + (2.0_dp*eps_x)*N_Ge; q_nsi_etau = (eps_y)*Z_Ge + (2.0_dp*eps_y)*N_Ge
        case(10);q_nsi_ee = (2.0_dp*eps_x)*Z_Ge + (eps_x)*N_Ge; q_nsi_etau = (eps_y)*Z_Ge + (2.0_dp*eps_y)*N_Ge
        case(11);q_nsi_ee = (eps_x)*Z_Ge + (2.0_dp*eps_x)*N_Ge; q_nsi_etau = (2.0_dp*eps_y)*Z_Ge + (eps_y)*N_Ge
        case(12);q_nsi_emu = (2.0_dp*eps_x)*Z_Ge + (eps_x)*N_Ge; q_nsi_etau = (2.0_dp*eps_y)*Z_Ge + (eps_y)*N_Ge
        case(13);q_nsi_emu = (eps_x)*Z_Ge + (2.0_dp*eps_x)*N_Ge; q_nsi_etau = (eps_y)*Z_Ge + (2.0_dp*eps_y)*N_Ge
        case(14);q_nsi_emu = (2.0_dp*eps_x)*Z_Ge + (eps_x)*N_Ge; q_nsi_etau = (eps_y)*Z_Ge + (2.0_dp*eps_y)*N_Ge
        case(15);q_nsi_emu = (eps_x)*Z_Ge + (2.0_dp*eps_x)*N_Ge; q_nsi_etau = (2.0_dp*eps_y)*Z_Ge + (eps_y)*N_Ge
        end select

        q_eff2 = (QW_SM + q_nsi_ee)**2 + q_nsi_emu**2 + q_nsi_etau**2
        R_th(:) = R_unit(:) * q_eff2

        ! Minimización analítica de alpha[cite: 1]
        sum1 = sum((R_exp * R_th) / sigma_exp**2); sum2 = sum(R_th**2 / sigma_exp**2)
        alpha_best = (sum1 - sum2) / (sum2 + 1.0_dp/sigma_alpha**2)
        chi2 = sum(((R_exp - (1.0_dp + alpha_best)*R_th)**2) / sigma_exp**2) + (alpha_best/sigma_alpha)**2
        
        write(u_out, '(3ES15.6)') eps_x, eps_y, chi2
     end do
  end do

  close(u_out)
  write(*,*) 'Archivo CONUS+ generado para ipar=', ipar

contains

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

end program chi2_nsi_2D
