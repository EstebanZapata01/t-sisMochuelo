!=======================================================================
! Programa: chi2_nsi_2D (Versión CONUS+ / Germanio)
!   - Usa matriz de resolución gaussiana para Germanio[cite: 6].
!   - Datos experimentales de excesos de CONUS+[cite: 2].
!   - Eje X (eps_x): Primera columna, usualmente sabor electrónico (ee).
!   - Eje Y (eps_y): Segunda columna, sabores mu o tau.
!
! Rol en la tesis: REFERENCIA validada. Reproduce el ajuste NSI 2D de
!   V. De Romeri, D. K. Papoulias, G. Sanchez Garcia, Phys. Rev. D 111,
!   075025 (2025) (Fig. 4). El analisis NSI de RED-100
!   (chi2red100_nest.f90 / chi2_ideal_nest.f90) es este mismo metodo
!   portado a Xe/Ar: mismos 15 casos ipar, misma carga efectiva Ec. (9),
!   misma minimizacion analitica de alpha. Lo que cambia: aqui hay DATOS
!   (senal a 3.7 sigma) -> region PERMITIDA; en RED-100 no -> sensibilidad.
! Tesis   : metodologia.tex Sec. 6 (CONUS+ frente a RED-100).
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
  ! ipar: 5 por defecto (sin argumentos, igual que siempre); opcionalmente
  ! por linea de comandos (ej. para los otros planos de la Fig. 4 de
  ! CONUS+). Con ipar!=5 el .dat/.txt de la grilla llevan el sufijo
  ! _ipar<N> para no pisar los de ipar=5 que usa el resto del pipeline;
  ! generic_input_conus.dat SI se sobreescribe siempre (es aditivo por
  ! diseno) -- volver a correr con ipar=5 al final restaura la validacion
  ! cruzada Xe<->CONUS+.
  block
    character(len=20) :: arg
    integer :: nargs
    ipar = 5
    nargs = command_argument_count()
    if (nargs >= 1) then
       call get_command_argument(1, arg)
       read(arg, *) ipar
    end if
  end block

  eps_y_min = -1.0_dp; eps_y_max = 1.0_dp
  eps_x_min = -1.0_dp; eps_x_max = 1.0_dp
  deps_y = (eps_y_max - eps_y_min) / (n_y - 1)
  deps_x = (eps_x_max - eps_x_min) / (n_x - 1)

  outdir = '../../datos/'
  if (ipar == 5) then
     filename = trim(outdir) // 'chi2_nsi_2Dconus.dat'
     configfile = trim(outdir) // 'nsi_config_conus.txt'
  else
     write(filename, '(A,A,I0,A)') trim(outdir), 'chi2_nsi_2Dconus_ipar', ipar, '.dat'
     write(configfile, '(A,A,I0,A)') trim(outdir), 'nsi_config_conus_ipar', ipar, '.txt'
  end if

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

  ! ================== 1b. VOLCADO ADITIVO para chi2_nsi_generic ==================
  ! Validacion cruzada Xe <-> CONUS+ con el MISMO motor chi2+NSI (ver
  ! metodologia.tex). NO altera ningun calculo/archivo anterior: solo agrega
  ! este archivo nuevo con los datos crudos de entrada. R_pred se guarda como
  ! R_unit(bin)*QW_SM^2, es decir la prediccion SM completa (eps=0), en la
  ! misma convencion que R_SM en chi2.f90 (Xe).
  block
    integer :: u_gen, ib
    open(newunit=u_gen, file=trim(outdir)//'generic_input_conus.dat', status='replace')
    write(u_gen,'(A)') '# Z  N  sigma_alpha  n_bins  ipar'
    write(u_gen,'(F8.2,F8.2,F8.4,2I6)') real(Z_Ge,dp), N_Ge, sigma_alpha, nbins, ipar
    write(u_gen,'(A)') '# bin  R_exp  sigma_exp  R_pred(=R_unit*QW_SM^2)'
    do ib = 1, nbins
       write(u_gen,'(I5,3ES16.7)') ib, R_exp(ib), sigma_exp(ib), R_unit(ib)*QW_SM**2
    end do
    close(u_gen)
    write(*,*) 'Volcado generic_input_conus.dat (validacion cruzada, aditivo)'
  end block

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
        write(u_out, '(3ES15.6)') eps_x, eps_y, chi2_en(eps_x, eps_y)
     end do
  end do

  close(u_out)
  write(*,*) 'Archivo CONUS+ generado para ipar=', ipar

  ! ================== 4. CORTES 1D FINOS, el otro parametro EXACTO en 0 ==================
  ! El corte 1D de la seccion 3 toma la fila/columna del grid 2D mas CERCANA
  ! a 0 (resolucion deps~0.002), no el otro parametro EXACTAMENTE en 0; eso
  ! da intervalos levemente asimetricos. Aqui se barre cada eje por separado
  ! con el otro en 0.0 exacto, igual que la Tabla II de De Romeri et al.
  block
    integer, parameter :: n_1d = 4000
    integer :: k
    real(dp) :: e, deps_1d
    character(len=250) :: f_1dx, f_1dy
    integer :: u1
    deps_1d = 2.0_dp / real(n_1d - 1, dp)
    if (ipar == 5) then
       f_1dx = trim(outdir) // 'chi2_nsi_1Dconus_x.dat'
       f_1dy = trim(outdir) // 'chi2_nsi_1Dconus_y.dat'
    else
       write(f_1dx, '(A,A,I0,A)') trim(outdir), 'chi2_nsi_1Dconus_ipar', ipar, '_x.dat'
       write(f_1dy, '(A,A,I0,A)') trim(outdir), 'chi2_nsi_1Dconus_ipar', ipar, '_y.dat'
    end if
    open(newunit=u1, file=trim(f_1dx), status='replace')
    write(u1, '(A)') '# eps_x (eps_y=0 exacto)   chi2'
    do k = 0, n_1d-1
       e = -1.0_dp + k*deps_1d
       write(u1, '(2ES15.6)') e, chi2_en(e, 0.0_dp)
    end do
    close(u1)
    open(newunit=u1, file=trim(f_1dy), status='replace')
    write(u1, '(A)') '# eps_y (eps_x=0 exacto)   chi2'
    do k = 0, n_1d-1
       e = -1.0_dp + k*deps_1d
       write(u1, '(2ES15.6)') e, chi2_en(0.0_dp, e)
    end do
    close(u1)
    write(*,*) 'Cortes 1D finos (eps_y=0 y eps_x=0 exactos) escritos.'
  end block

contains

  function chi2_en(ex, ey) result(chi2_out)
    real(dp), intent(in) :: ex, ey
    real(dp) :: chi2_out, qee, qemu, qetau, qeff2, Rth(nbins), s1, s2, abest
    qee = 0.0_dp; qemu = 0.0_dp; qetau = 0.0_dp
    select case(ipar)
    case(1); qee = (2.0_dp*ey + ex)*Z_Ge + (ey + 2.0_dp*ex)*N_Ge
    case(2); qemu = (2.0_dp*ey + ex)*Z_Ge + (ey + 2.0_dp*ex)*N_Ge
    case(3); qetau = (2.0_dp*ey + ex)*Z_Ge + (ey + 2.0_dp*ex)*N_Ge
    case(4); qee = (2.0_dp*ex)*Z_Ge + (ex)*N_Ge; qemu = (2.0_dp*ey)*Z_Ge + (ey)*N_Ge
    case(5); qee = (ex)*Z_Ge + (2.0_dp*ex)*N_Ge; qemu = (ey)*Z_Ge + (2.0_dp*ey)*N_Ge
    case(6); qee = (2.0_dp*ex)*Z_Ge + (ex)*N_Ge; qemu = (ey)*Z_Ge + (2.0_dp*ey)*N_Ge
    case(7); qee = (ex)*Z_Ge + (2.0_dp*ex)*N_Ge; qemu = (2.0_dp*ey)*Z_Ge + (ey)*N_Ge
    case(8); qee = (2.0_dp*ex)*Z_Ge + (ex)*N_Ge; qetau = (2.0_dp*ey)*Z_Ge + (ey)*N_Ge
    case(9); qee = (ex)*Z_Ge + (2.0_dp*ex)*N_Ge; qetau = (ey)*Z_Ge + (2.0_dp*ey)*N_Ge
    case(10);qee = (2.0_dp*ex)*Z_Ge + (ex)*N_Ge; qetau = (ey)*Z_Ge + (2.0_dp*ey)*N_Ge
    case(11);qee = (ex)*Z_Ge + (2.0_dp*ex)*N_Ge; qetau = (2.0_dp*ey)*Z_Ge + (ey)*N_Ge
    case(12);qemu = (2.0_dp*ex)*Z_Ge + (ex)*N_Ge; qetau = (2.0_dp*ey)*Z_Ge + (ey)*N_Ge
    case(13);qemu = (ex)*Z_Ge + (2.0_dp*ex)*N_Ge; qetau = (ey)*Z_Ge + (2.0_dp*ey)*N_Ge
    case(14);qemu = (2.0_dp*ex)*Z_Ge + (ex)*N_Ge; qetau = (ey)*Z_Ge + (2.0_dp*ey)*N_Ge
    case(15);qemu = (ex)*Z_Ge + (2.0_dp*ex)*N_Ge; qetau = (2.0_dp*ey)*Z_Ge + (ey)*N_Ge
    end select
    qeff2 = (QW_SM + qee)**2 + qemu**2 + qetau**2
    Rth = R_unit * qeff2
    s1 = sum((R_exp * Rth) / sigma_exp**2); s2 = sum(Rth**2 / sigma_exp**2)
    abest = (s1 - s2) / (s2 + 1.0_dp/sigma_alpha**2)
    chi2_out = sum(((R_exp - (1.0_dp + abest)*Rth)**2) / sigma_exp**2) + (abest/sigma_alpha)**2
  end function chi2_en

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
