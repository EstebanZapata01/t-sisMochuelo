!=======================================================================
! Programa: chi2red100  (matriz NSI 2D - RAMA UNIFICADA con chi2.f90)
!
!   La NSI entra en CEvNS solo a traves de la carga efectiva:
!       q_eff2(eps) = (Q_w + q_ee)^2 + q_emu^2 + q_etau^2
!   que reescala TODO el espectro igual que la amplitud A de chi2.f90.
!   => el mapa NSI 2D es chi2.f90 evaluado en  A = q_eff2(eps)/Q_w^2
!      en cada punto de la grilla, reutilizando el ajuste ON-OFF real
!      (datos digitalizados + errores estadisticos reales, que ya
!      reproducen la sensibilidad ~46xSM del paper arXiv:2411.18641).
!
!   chi2(A, alpha) = Sum_i (dNi - A*(1+alpha)*Ri)^2/si^2 + (alpha/sigF)^2
!     minimizado analiticamente sobre alpha (nuisance de flujo 16.9%).
!
!   Entrada : ionization_spectra_detallado.dat  (espectro PE teorico, red100PE.f90)
!   Salida  : chi2_nsi_2DXe.dat  (eps_x  eps_y  chi2)
!             nsi_configXe.txt   (etiquetas de los ejes)
!
!   Caso por defecto ipar=5: eps_ee^dV (x) vs eps_emu^dV (y)
!
! Pipeline: rama NSI 2D CON DATOS REALES (residuo ON-OFF 2024). Es la
!           version "no ideal" de chi2_ideal_nest.f90.
! Tesis   : metodologia.tex Sec. 5 (formalismo NSI) y Sec. 6 (CONUS+ vs
!           RED-100). El mapeo eps -> q_nsi es IDENTICO a 2pchi2.f90
!           (CONUS+); los 15 casos ipar deben conservarse.
! Decision metodologica clave: A_amp(eps) = q_eff2/Q_w^2 reescala todos
!   los bins por igual (las NSI vectoriales no distorsionan el espectro).
!=======================================================================
program chi2red100
  use constants, only: dp, Z_Ge, N_Ge, eff_ROI
  implicit none

  ! ----- dataset ON-OFF digitalizado (MISMO que chi2.f90; mantener en sync) -----
  integer, parameter :: n_datos = 15
  real(dp) :: pe_dat(n_datos), dN_dat(n_datos), sigma_dat(n_datos)

  ! ----- prediccion SM (espectro PE teorico de red100PE.f90) -----
  integer :: n_pred
  real(dp), allocatable :: pe_pred(:), R_pred(:)
  real(dp) :: R_SM(n_datos)

  ! ----- analisis -----
  real(dp), parameter :: sin2th_SM  = 0.23857_dp
  real(dp), parameter :: sigma_F    = 0.169_dp        ! nuisance de normalizacion de flujo
  integer,  parameter :: ipar       = 5              ! <-- caso NSI
  integer,  parameter :: n_u = 1000, n_d = 1000

  real(dp) :: S1, S2, S3, QW_SM
  real(dp) :: q_nsi_ee, q_nsi_emu, q_nsi_etau, q_eff2, A_amp, alpha_best, chi2
  real(dp) :: eps_x, eps_y, eps_min, eps_max, deps
  real(dp) :: slope, cols(9)
  integer  :: i, j, j0, ios, u_pred, u_out, u_conf
  ! --- validacion ---
  real(dp) :: chi2_A1, chi2_1d_min, A_1d_min, A90_eff, Aeff_lo, Aeff_hi, Aa, cc
  integer  :: n_in90, n_in68, n_tot
  character(len=256) :: line
  character(len=250) :: datadir, f_pred, f_out, f_conf
  character(len=60)  :: xlabel, ylabel

  datadir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  f_pred  = trim(datadir)//'ionization_spectra_detallado.dat'
  f_out   = trim(datadir)//'chi2_nsi_2DXe.dat'
  f_conf  = trim(datadir)//'nsi_configXe.txt'

  ! ==================================================================
  ! 1. DATOS DIGITALIZADOS DEL ARTICULO (residuo ON-OFF, counts/kg/dia)
  ! ==================================================================
  pe_dat(1)  = 112.60128_dp; dN_dat(1)  = -0.21102_dp; sigma_dat(1)  = 0.20630_dp
  pe_dat(2)  = 117.86781_dp; dN_dat(2)  =  0.09291_dp; sigma_dat(2)  = 0.18583_dp
  pe_dat(3)  = 123.14301_dp; dN_dat(3)  = -0.29449_dp; sigma_dat(3)  = 0.18583_dp
  pe_dat(4)  = 128.34364_dp; dN_dat(4)  =  0.10394_dp; sigma_dat(4)  = 0.14961_dp
  pe_dat(5)  = 133.57848_dp; dN_dat(5)  =  0.14173_dp; sigma_dat(5)  = 0.13071_dp
  pe_dat(6)  = 138.88895_dp; dN_dat(6)  = -0.01732_dp; sigma_dat(6)  = 0.12756_dp
  pe_dat(7)  = 144.19912_dp; dN_dat(7)  = -0.17323_dp; sigma_dat(7)  = 0.12913_dp
  pe_dat(8)  = 149.46939_dp; dN_dat(8)  =  0.09134_dp; sigma_dat(8)  = 0.10551_dp
  pe_dat(9)  = 154.71783_dp; dN_dat(9)  = -0.01417_dp; sigma_dat(9)  = 0.09291_dp
  pe_dat(10) = 159.96419_dp; dN_dat(10) = -0.09764_dp; sigma_dat(10) = 0.09291_dp
  pe_dat(11) = 165.24581_dp; dN_dat(11) =  0.04724_dp; sigma_dat(11) = 0.06929_dp
  pe_dat(12) = 170.54956_dp; dN_dat(12) = -0.04094_dp; sigma_dat(12) = 0.07402_dp
  pe_dat(13) = 175.78321_dp; dN_dat(13) =  0.00945_dp; sigma_dat(13) = 0.06614_dp
  pe_dat(14) = 181.02373_dp; dN_dat(14) = -0.01260_dp; sigma_dat(14) = 0.06772_dp
  pe_dat(15) = 186.30357_dp; dN_dat(15) =  0.15118_dp; sigma_dat(15) = 0.05512_dp

  ! ==================================================================
  ! 2. LEER PREDICCION SM (9 col: PE Total 1SE..7SE) y pesar la ROI
  !    R_pred(PE) = Sum_{k=4}^{7} col_kSE(PE) * eff_ROI(k)
  ! ==================================================================
  n_pred = 0
  open(newunit=u_pred, file=f_pred, status='old', action='read')
  do
    read(u_pred, '(A)', iostat=ios) line
    if (ios < 0) exit
    if (ios > 0) cycle
    line = adjustl(line)
    if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
    n_pred = n_pred + 1
  end do
  rewind(u_pred)
  allocate(pe_pred(n_pred), R_pred(n_pred))

  j = 0
  do
    read(u_pred, '(A)', iostat=ios) line
    if (ios < 0) exit
    if (ios > 0) cycle
    line = adjustl(line)
    if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
    read(line, *, iostat=ios) cols
    if (ios /= 0) cycle
    j = j + 1
    pe_pred(j) = cols(1)
    R_pred(j)  = cols(6)*eff_ROI(4) + cols(7)*eff_ROI(5) &
               + cols(8)*eff_ROI(6) + cols(9)*eff_ROI(7)
  end do
  close(u_pred)
  write(*,'(A,I5,A)') '  Leidos ', n_pred, ' bins PE de la prediccion SM'

  ! ==================================================================
  ! 3. INTERPOLAR PREDICCION A LOS BINS DE DATOS
  ! ==================================================================
  do i = 1, n_datos
    R_SM(i) = 0.0_dp
    j0 = -1
    do j = 1, n_pred - 1
      if (pe_pred(j) <= pe_dat(i) .and. pe_dat(i) <= pe_pred(j+1)) then
        j0 = j; exit
      end if
    end do
    if (j0 > 0) then
      slope   = (R_pred(j0+1) - R_pred(j0)) / (pe_pred(j0+1) - pe_pred(j0))
      R_SM(i) = R_pred(j0) + slope * (pe_dat(i) - pe_pred(j0))
    else
      write(*,'(A,F7.1)') '  AVISO: bin PE fuera de rango -> R_SM=0 en PE=', pe_dat(i)
    end if
  end do

  ! ==================================================================
  ! 4. SUMAS AUXILIARES (independientes de A)
  ! ==================================================================
  S1 = 0.0_dp; S2 = 0.0_dp; S3 = 0.0_dp
  do i = 1, n_datos
    S1 = S1 + dN_dat(i) * R_SM(i)  / sigma_dat(i)**2
    S2 = S2 + R_SM(i)**2           / sigma_dat(i)**2
    S3 = S3 + dN_dat(i)**2         / sigma_dat(i)**2
  end do

  ! Carga debil SM (convencion con prefactor 1/pi de xsections_nest)
  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*sin2th_SM)/2.0_dp * Z_Ge

  ! ==================================================================
  ! 4b. VALIDACION PRE-BARRIDO (perfil 1D en la amplitud A)
  ! ==================================================================
  chi2_A1 = chi2_of_A(1.0_dp)                 ! chi2 en el punto SM (A=1)

  chi2_1d_min = 1.0e30_dp; A_1d_min = 0.0_dp
  Aa = 0.0_dp
  do
     cc = chi2_of_A(Aa)
     if (cc < chi2_1d_min) then; chi2_1d_min = cc; A_1d_min = Aa; end if
     Aa = Aa + 0.02_dp
     if (Aa > 400.0_dp) exit
  end do
  ! cruce Delta chi2 = 2.706 (90% 1 g.d.l.) subiendo desde A_1d_min
  A90_eff = -1.0_dp; Aa = A_1d_min
  do
     Aa = Aa + 0.02_dp
     if (Aa > 400.0_dp) exit
     if (chi2_of_A(Aa) - chi2_1d_min >= 2.706_dp) then; A90_eff = Aa; exit; end if
  end do

  write(*,'(/,A)') '=============== VALIDACION chi2red100_nest (rama ON-OFF, vs paper) ==============='
  write(*,'(A,F12.5)')   ' Q_w                              = ', QW_SM
  write(*,'(A,F12.3)')   ' Q_w^2                            = ', QW_SM**2
  write(*,'(A,I3)')      ' Bins de datos ON-OFF             = ', n_datos
  write(*,'(A,ES12.4)')  ' Sum R_SM (prediccion en la ROI)  = ', sum(R_SM(1:n_datos))
  write(*,'(A,3ES12.4)') ' S1, S2, S3                       = ', S1, S2, S3
  write(*,'(A,F10.4)')   ' chi2_min (perfil 1D en A)        = ', chi2_1d_min
  write(*,'(A,F10.4,A)') ' chi2_min / ndof                  = ', chi2_1d_min/real(n_datos,dp), &
       '   (ndof = 15)'
  write(*,'(A,F10.4)')   ' A en el minimo (x SM)            = ', A_1d_min
  write(*,'(A,F10.4)')   ' Delta chi2 en el SM (A=1)        = ', chi2_A1 - chi2_1d_min
  write(*,'(A,F10.3,A)') ' A_90 efectivo (Dchi2=2.706)      = ', A90_eff, &
       '  x SM   (comparar chi2.f90 ~111; paper KI 94, SM2018 63)'

  ! ==================================================================
  ! 5. ETIQUETAS DE LOS EJES SEGUN ipar
  ! ==================================================================
  select case(ipar)
  case(1);  xlabel='$\epsilon_{ee}^{dV}$';   ylabel='$\epsilon_{ee}^{uV}$'
  case(2);  xlabel='$\epsilon_{e\mu}^{dV}$';  ylabel='$\epsilon_{e\mu}^{uV}$'
  case(3);  xlabel='$\epsilon_{e\tau}^{dV}$'; ylabel='$\epsilon_{e\tau}^{uV}$'
  case(4);  xlabel='$\epsilon_{ee}^{uV}$';   ylabel='$\epsilon_{e\mu}^{uV}$'
  case(5);  xlabel='$\epsilon_{ee}^{dV}$';   ylabel='$\epsilon_{e\mu}^{dV}$'
  case(6);  xlabel='$\epsilon_{ee}^{uV}$';   ylabel='$\epsilon_{e\mu}^{dV}$'
  case(7);  xlabel='$\epsilon_{ee}^{dV}$';   ylabel='$\epsilon_{e\mu}^{uV}$'
  case(8);  xlabel='$\epsilon_{ee}^{uV}$';   ylabel='$\epsilon_{e\tau}^{uV}$'
  case(9);  xlabel='$\epsilon_{ee}^{dV}$';   ylabel='$\epsilon_{e\tau}^{dV}$'
  case(10); xlabel='$\epsilon_{ee}^{uV}$';   ylabel='$\epsilon_{e\tau}^{dV}$'
  case(11); xlabel='$\epsilon_{ee}^{dV}$';   ylabel='$\epsilon_{e\tau}^{uV}$'
  case(12); xlabel='$\epsilon_{e\mu}^{uV}$';  ylabel='$\epsilon_{e\tau}^{uV}$'
  case(13); xlabel='$\epsilon_{e\mu}^{dV}$';  ylabel='$\epsilon_{e\tau}^{dV}$'
  case(14); xlabel='$\epsilon_{e\mu}^{uV}$';  ylabel='$\epsilon_{e\tau}^{dV}$'
  case(15); xlabel='$\epsilon_{e\mu}^{dV}$';  ylabel='$\epsilon_{e\tau}^{uV}$'
  end select
  open(newunit=u_conf, file=f_conf, status='replace')
  write(u_conf, '(A)') trim(xlabel); write(u_conf, '(A)') trim(ylabel)
  close(u_conf)

  ! ==================================================================
  ! 6. BARRIDO NSI 2D
  ! ==================================================================
  eps_min = -1.0_dp; eps_max = 1.0_dp
  deps = (eps_max - eps_min) / real(n_u - 1, dp)

  open(newunit=u_out, file=f_out, status='replace')
  write(u_out, '(A)') '# eps_x   eps_y   chi2'
  write(*,*) 'Iniciando barrido NSI (1000x1000) - rama unificada ON-OFF...'

  n_in90 = 0; n_in68 = 0; n_tot = 0
  Aeff_lo = 1.0e30_dp; Aeff_hi = -1.0e30_dp

  do i = 0, n_u-1
     eps_y = eps_min + i * deps
     do j = 0, n_d-1
        eps_x = eps_min + j * deps
        q_nsi_ee = 0.0_dp; q_nsi_emu = 0.0_dp; q_nsi_etau = 0.0_dp

        select case(ipar)
        case(1); q_nsi_ee   = (2.0_dp*eps_y + eps_x)*Z_Ge + (eps_y + 2.0_dp*eps_x)*N_Ge
        case(2); q_nsi_emu  = (2.0_dp*eps_y + eps_x)*Z_Ge + (eps_y + 2.0_dp*eps_x)*N_Ge
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
        A_amp  = q_eff2 / QW_SM**2          ! amplitud CEvNS efectiva (A=1 -> SM)

        ! minimo analitico sobre el nuisance de flujo alpha
        alpha_best = A_amp*(S1 - A_amp*S2) / (1.0_dp/sigma_F**2 + A_amp**2*S2)

        chi2 = 0.0_dp
        do j0 = 1, n_datos
           chi2 = chi2 + ((dN_dat(j0) - A_amp*(1.0_dp + alpha_best)*R_SM(j0)) / sigma_dat(j0))**2
        end do
        chi2 = chi2 + (alpha_best/sigma_F)**2

        n_tot = n_tot + 1
        if (chi2 - chi2_1d_min <= 4.605_dp) n_in90 = n_in90 + 1   ! 90% 2 g.d.l.
        if (chi2 - chi2_1d_min <= 2.300_dp) n_in68 = n_in68 + 1   ! 68% 2 g.d.l.
        if (A_amp < Aeff_lo) Aeff_lo = A_amp
        if (A_amp > Aeff_hi) Aeff_hi = A_amp

        write(u_out, '(3ES15.6)') eps_x, eps_y, chi2
     end do
     write(u_out, *)
  end do
  close(u_out)

  write(*,'(/,A)') ' --- Barrido NSI (ipar=5: eps_ee^dV vs eps_emu^dV, |eps| <= 1) ---'
  write(*,'(A,2F9.3)') ' Amplitud efectiva q_eff2/Q_w^2 en la grilla: min, max = ', Aeff_lo, Aeff_hi
  write(*,'(A,F8.3,A)') ' Fraccion del plano con Dchi2 < 4.605 (90%, 2 gdl) = ', &
       real(n_in90,dp)/real(n_tot,dp), '   (1.0 => sin restriccion)'
  write(*,'(A,F8.3)')   ' Fraccion del plano con Dchi2 < 2.300 (68%, 2 gdl) = ', &
       real(n_in68,dp)/real(n_tot,dp)
  write(*,'(A,I0,A)')  ' --- Matriz NSI generada (ipar=', ipar, ') - rama ON-OFF unificada ---'
  write(*,'(A,/)') '================================================================================'

  deallocate(pe_pred, R_pred)

contains

  ! chi2 del ajuste ON-OFF para una amplitud CEvNS A (min analitico sobre alpha)
  function chi2_of_A(A) result(c)
    real(dp), intent(in) :: A
    real(dp) :: c, al
    integer  :: m
    al = A*(S1 - A*S2) / (1.0_dp/sigma_F**2 + A**2*S2)
    c = 0.0_dp
    do m = 1, n_datos
       c = c + ((dN_dat(m) - A*(1.0_dp + al)*R_SM(m)) / sigma_dat(m))**2
    end do
    c = c + (al/sigma_F)**2
  end function chi2_of_A

end program chi2red100
