!=======================================================================
! Programa: chi2_ideal_nest   (COMPARACION IDEAL SIMETRICA Xe vs Ar)
!
!   MISMO codigo para los dos blancos (solo cambia el parametro TAG y,
!   por la carpeta, el modulo `constants` y la tabla NEST que lee
!   `mod_tnr_to_e`). Sirve para una comparacion Xe-vs-Ar limpia:
!
!     - Sensibilidad ASIMOV DE CONTEO PURO en el ROI (sin fondo, sin
!       sistematicos, sin eficiencia de seleccion: eff_ROI = 1):
!           dN_k = R_k ,  sigma_k = sqrt(N_k)
!           chi2(A) = (1 - A)^2 * N_tot_ROI ,  A_best = 1
!           A_90 (1 gdl) : A_90 = 1 + sqrt(2.706 / N_tot_ROI)
!       Mide la RESPUESTA INTRINSECA del blanco
!       (sigma_CEvNS . N^2 . charge yield . umbral). NO es la
!       sensibilidad experimental alcanzable.
!
!     - VENTANA SIMETRICA: se barre el umbral inferior N_e >= {1,2,3,4}
!       (el corte N_e>=4 de RED-100 es fondo de electron unico, NO
!       sensibilidad de los aparatos) y NO se pone corte superior
!       (N_e <= NE_HI = 30; la cola dura del Ar llega a ~N_e 30-40 y no
!       hay que tirarla). Cualquier corte comun es legitimo mientras sea
!       el mismo para los dos blancos.
!
!     - FLUCTUACION F(T): se corre con la F de NEST propia del blanco
!       (de la tabla, columna 3) y con F = 1 (Poisson) como sistematico.
!
!     - NSI 2D: la NSI entra solo por
!           q_eff2(eps) = (Q_W + q_ee)^2 + q_emu^2 + q_etau^2
!           A_amp(eps)  = q_eff2 / Q_W^2
!           Delta chi2  = (1 - A_amp)^2 * N_tot_ROI
!       El par de parametros lo fija `ipar` (1..15, mismo mapeo que
!       chi2red100_nest.f90 / 2pchi2.f90 de CONUS+); por defecto ipar=5
!       (eps_ee^dV vs eps_emu^dV). Se hace para N_e>=1 y para N_e>=4
!       (F de NEST). Las etiquetas de los ejes se escriben en
!       nsi_config_ideal_<TAG>.txt.
!
!   El limite OBSERVADO real de Xe (~111 xSM, rama chi2.f90 con datos
!   ON-OFF 2024 y eff_ROI digitalizada) se reporta APARTE.
!   Para Ar solo es fisico asumiendo ARGON DEPLETADO (UAr).
!
!   Entrada : nada (calcula el espectro SM internamente)
!             lee la tabla NEST via mod_tnr_to_e
!   Salidas (datos/):
!     espectro_Ne_ideal_<TAG>.dat  : N_e  R_bin_Fnest  R_bin_F1   [ev/(kg dia)]
!     sensib_ideal_<TAG>.dat       : NE_LO F_mode R_tot frac Ntot@192 A_90(x1..x335)
!     chi2_nsi_2D<TAG>_ideal.dat       : eps_x eps_y Dchi2   (N_e>=1, F de NEST)
!     chi2_nsi_2D<TAG>_ideal_ne4.dat   : eps_x eps_y Dchi2   (N_e>=4, F de NEST)
!     nsi_config_ideal_<TAG>.txt   : etiqueta LaTeX del eje x / eje y (segun ipar)
!=======================================================================
program chi2_ideal_nest
  use constants
  use mod_tnr_to_e
  use xsections_nest
  use mod_stats,   only: binomial_prob
  use flux,        only: flujo_diferencial, E_nu_max
  implicit none

  character(len=*), parameter :: TAG = 'Ar'   ! <-- unica diferencia con la copia de Ar

  ! ----- grilla de integracion -----
  integer,  parameter :: n_T = 800, n_E = 2000
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, peso
  real(dp) :: tasa_Comb, QW_SM, lambda
  integer  :: i_T, i_E, k, n_F, im
  real(dp) :: p_F

  ! ----- espectro N_e (ROI sin corte superior) -----
  integer,  parameter :: NE_HI = 30
  real(dp) :: array_tasa_Comb(n_T)
  real(dp) :: R_bin(NE_HI, 2)          ! (N_e, F_mode)  tasa SM  [ev/(kg dia)]
  real(dp) :: R_tot, frac, Ntot_roi, A90
  integer,  parameter :: NE_LO_scan(4) = (/ 1, 2, 3, 4 /)
  integer  :: is
  character(len=8) :: fname(2) = (/ 'Fnest   ', 'F1      ' /)

  ! ----- Ntot_roi guardado para el barrido NSI (F de NEST, im=1) -----
  real(dp) :: Ntot_ne1, Ntot_ne4

  ! ----- NSI 2D -----
  integer,  parameter :: ipar = 5     ! <-- par de parametros NSI (1..15); ver select case
  integer,  parameter :: n_u = 1000, n_d = 1000
  real(dp) :: eps_x, eps_y, eps_min, eps_max, deps
  real(dp) :: q_nsi_ee, q_nsi_emu, q_nsi_etau, q_eff2, A_amp, dchi2
  real(dp) :: Aeff_lo, Aeff_hi
  integer  :: i, j, u, n_in90, n_in68, n_tot
  character(len=40) :: xlabel, ylabel

  ! ----- exposicion -----
  real(dp), parameter :: mult(7) = &
       (/ 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 50.0_dp, 100.0_dp, 335.0_dp /)

  character(len=300) :: datadir, f_esp, f_sens, f_nsi1, f_nsi4, f_conf
  real(dp) :: dummy

  datadir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  f_esp  = trim(datadir)//'espectro_Ne_ideal_'//TAG//'.dat'
  f_sens = trim(datadir)//'sensib_ideal_'//TAG//'.dat'
  f_nsi1 = trim(datadir)//'chi2_nsi_2D'//TAG//'_ideal.dat'
  f_nsi4 = trim(datadir)//'chi2_nsi_2D'//TAG//'_ideal_ne4.dat'
  f_conf = trim(datadir)//'nsi_config_ideal_'//TAG//'.txt'

  ! ----- etiquetas de los ejes segun ipar (mismo mapeo que chi2red100_nest.f90) -----
  select case(ipar)
  case(1);  xlabel='$\epsilon_{ee}^{dV}$';     ylabel='$\epsilon_{ee}^{uV}$'
  case(2);  xlabel='$\epsilon_{e\mu}^{dV}$';    ylabel='$\epsilon_{e\mu}^{uV}$'
  case(3);  xlabel='$\epsilon_{e\tau}^{dV}$';   ylabel='$\epsilon_{e\tau}^{uV}$'
  case(4);  xlabel='$\epsilon_{ee}^{uV}$';      ylabel='$\epsilon_{e\mu}^{uV}$'
  case(5);  xlabel='$\epsilon_{ee}^{dV}$';      ylabel='$\epsilon_{e\mu}^{dV}$'
  case(6);  xlabel='$\epsilon_{ee}^{uV}$';      ylabel='$\epsilon_{e\mu}^{dV}$'
  case(7);  xlabel='$\epsilon_{ee}^{dV}$';      ylabel='$\epsilon_{e\mu}^{uV}$'
  case(8);  xlabel='$\epsilon_{ee}^{uV}$';      ylabel='$\epsilon_{e\tau}^{uV}$'
  case(9);  xlabel='$\epsilon_{ee}^{dV}$';      ylabel='$\epsilon_{e\tau}^{dV}$'
  case(10); xlabel='$\epsilon_{ee}^{uV}$';      ylabel='$\epsilon_{e\tau}^{dV}$'
  case(11); xlabel='$\epsilon_{ee}^{dV}$';      ylabel='$\epsilon_{e\tau}^{uV}$'
  case(12); xlabel='$\epsilon_{e\mu}^{uV}$';    ylabel='$\epsilon_{e\tau}^{uV}$'
  case(13); xlabel='$\epsilon_{e\mu}^{dV}$';    ylabel='$\epsilon_{e\tau}^{dV}$'
  case(14); xlabel='$\epsilon_{e\mu}^{uV}$';    ylabel='$\epsilon_{e\tau}^{dV}$'
  case(15); xlabel='$\epsilon_{e\mu}^{dV}$';    ylabel='$\epsilon_{e\tau}^{uV}$'
  case default
     write(*,*) 'ERROR: ipar fuera de rango (1..15): ', ipar
     stop 1
  end select

  open(newunit=u, file=f_conf, status='replace')
  write(u,'(A)') trim(xlabel)
  write(u,'(A)') trim(ylabel)
  close(u)

  call inicializar_nest()
  dummy = flujo_diferencial(1.0_dp)

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2   = QW_SM**2

  ! ==================================================================
  ! 1. ESPECTRO DE RETROCESO CEvNS SM   [ev/(kg dia keV)]
  ! ==================================================================
  T_nr_min = 0.20_dp / 1000.0_dp
  T_nr_max = 3.5_dp  / 1000.0_dp     ! ref.[46] + cinematica E_nu=8 MeV (T_max=3.44 keV);
                                     ! trunca la cola E_nu 8-10 MeV (<1 %). Antes 6.0. En Xe el
                                     ! corte cinematico anula el integrando > ~1.6 keV
  dT = (T_nr_max - T_nr_min) / (n_T - 1); dT_keV = dT * 1000.0_dp
  dE = E_nu_max / (n_E - 1)

  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Comb = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle
        peso = merge(0.5_dp, 1.0_dp, i_E == 1 .or. i_E == n_E)
        tasa_Comb = tasa_Comb + flujo_diferencial(E_nu) * dsigma_dT(E_nu, T_nr) * dE * peso
     end do
     array_tasa_Comb(i_T) = tasa_Comb * (NA / A_Ge) * 1000.0_dp * 86400.0_dp
  end do

  ! ==================================================================
  ! 2. ESPECTRO EN N_e  (IDEAL: sin eff_ROI, sin livetime)
  !    im=1 -> F de NEST (tabla) ; im=2 -> F = 1 (Poisson)
  ! ==================================================================
  R_bin = 0.0_dp
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)

     do im = 1, 2
        if (im == 1) then
           call obtener_nest_binomial(T_nr * 1000.0_dp, n_F, p_F)
        else
           lambda = obtener_electrones_creados(T_nr * 1000.0_dp)   ! <N_e> creados
           p_F = 0.02_dp                                           ! limite Poisson
           n_F = nint(lambda / p_F)
           if (n_F < 1) n_F = 1
        end if
        do k = 1, NE_HI
           R_bin(k, im) = R_bin(k, im) + array_tasa_Comb(i_T) * dT_keV * peso * &
                          binomial_prob(k, n_F, p_F * EEE)
        end do
     end do
  end do

  ! ---- escribir el espectro N_e ----
  open(newunit=u, file=f_esp, status='replace')
  write(u,'(A)') '# N_e   R_bin_Fnest[ev/(kg dia)]   R_bin_F1[ev/(kg dia)]'
  do k = 1, NE_HI
     write(u,'(I4,2ES18.8)') k, R_bin(k,1), R_bin(k,2)
  end do
  close(u)

  ! ==================================================================
  ! 3. SENSIBILIDAD Asimov de conteo puro: barrido de umbral y de F
  ! ==================================================================
  write(*,'(/,A)') '======== COMPARACION IDEAL  ('//TAG//')  -- Asimov conteo puro ========'
  write(*,'(A,F12.5)')  ' Q_w                 = ', QW_SM
  write(*,'(A,F12.3)')  ' Q_w^2               = ', QW_SM**2
  write(*,'(A,F8.4)')   ' EEE                 = ', EEE
  write(*,'(A,F8.1,A)') ' Exposicion base     = ', exposure_ON_kgd, ' kg*dia'
  write(*,'(A,I0)')     ' NE_HI (corte sup.)  = ', NE_HI
  write(*,'(A,ES13.5)') ' R_tot 1..NE_HI (F nest) [ev/(kg dia)] = ', sum(R_bin(:,1))
  write(*,'(A,ES13.5)') ' R_tot 1..NE_HI (F=1)    [ev/(kg dia)] = ', sum(R_bin(:,2))

  open(newunit=u, file=f_sens, status='replace')
  write(u,'(A)') '# NE_LO  F_mode   R_tot_ROI[ev/kgd]   frac_vs_NE_LO1   Ntot_ROI@192   '// &
                 'A_90(x1)  A_90(x2)  A_90(x5)  A_90(x10)  A_90(x50)  A_90(x100)  A_90(x335)'

  do im = 1, 2
     write(*,'(/,A)') ' --- F_mode = '//trim(fname(im))//' ---'
     write(*,'(A)')   '   NE_LO   R_tot_ROI     frac      Ntot@192      A_90(base)'
     do is = 1, size(NE_LO_scan)
        R_tot = 0.0_dp
        do k = NE_LO_scan(is), NE_HI
           R_tot = R_tot + R_bin(k, im)
        end do
        frac     = R_tot / max(sum(R_bin(:,im)), 1.0e-300_dp)
        Ntot_roi = R_tot * exposure_ON_kgd
        A90      = 1.0_dp + sqrt(2.706_dp / max(Ntot_roi, 1.0e-300_dp))

        write(*,'(I6,ES14.5,F10.4,ES14.5,F12.4)') &
             NE_LO_scan(is), R_tot, frac, Ntot_roi, A90
        write(u,'(I6,3X,A8,ES18.6,F14.6,ES16.6)', advance='no') &
             NE_LO_scan(is), fname(im), R_tot, frac, Ntot_roi
        do i = 1, size(mult)
           write(u,'(F11.5)', advance='no') &
                1.0_dp + sqrt(2.706_dp / max(Ntot_roi*mult(i), 1.0e-300_dp))
        end do
        write(u,*)

        if (im == 1 .and. NE_LO_scan(is) == 1) Ntot_ne1 = Ntot_roi
        if (im == 1 .and. NE_LO_scan(is) == 4) Ntot_ne4 = Ntot_roi
     end do
  end do
  close(u)

  ! ==================================================================
  ! 4. BARRIDO NSI 2D  (ipar=5)  para N_e>=1 y N_e>=4  (F de NEST)
  ! ==================================================================
  eps_min = -1.0_dp; eps_max = 1.0_dp
  deps = (eps_max - eps_min) / real(n_u - 1, dp)

  call barrido_nsi(f_nsi1, Ntot_ne1, 1)
  call barrido_nsi(f_nsi4, Ntot_ne4, 4)

  write(*,'(/,A)') '================================================================================'

contains

  subroutine barrido_nsi(fname_out, Ntot, ne_lo_lab)
    character(len=*), intent(in) :: fname_out
    real(dp),         intent(in) :: Ntot
    integer,          intent(in) :: ne_lo_lab
    integer :: uu

    open(newunit=uu, file=fname_out, status='replace')
    write(uu,'(A,I0,A)') '# eps_x   eps_y   Delta_chi2   (IDEAL conteo puro, N_e>=', &
                         ne_lo_lab, ', F de NEST, exposicion base)'
    n_in90 = 0; n_in68 = 0; n_tot = 0
    Aeff_lo = 1.0e30_dp; Aeff_hi = -1.0e30_dp

    do i = 0, n_u-1
       eps_y = eps_min + i * deps
       do j = 0, n_d-1
          eps_x = eps_min + j * deps
          q_nsi_ee = 0.0_dp; q_nsi_emu = 0.0_dp; q_nsi_etau = 0.0_dp

          ! mapeo eps -> carga NSI, IDENTICO a chi2red100_nest.f90 / 2pchi2.f90
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
          A_amp  = q_eff2 / QW_SM**2
          dchi2  = (1.0_dp - A_amp)**2 * Ntot

          n_tot = n_tot + 1
          if (dchi2 <= 4.605_dp) n_in90 = n_in90 + 1
          if (dchi2 <= 2.300_dp) n_in68 = n_in68 + 1
          if (A_amp < Aeff_lo) Aeff_lo = A_amp
          if (A_amp > Aeff_hi) Aeff_hi = A_amp

          write(uu,'(3ES15.6)') eps_x, eps_y, dchi2
       end do
       write(uu,*)
    end do
    close(uu)

    write(*,'(/,A,I0,A,I0,A)') ' --- NSI 2D ('//TAG//', ipar=', ipar, &
         ': '//trim(xlabel)//' vs '//trim(ylabel)//', |eps|<=1), N_e>=', ne_lo_lab, ' ---'
    write(*,'(A,ES13.5)') '   N_tot_ROI usado            = ', Ntot
    write(*,'(A,2F10.3)') '   A_amp min, max en la caja  = ', Aeff_lo, Aeff_hi
    write(*,'(A,F8.4,A)') '   frac Dchi2 < 4.605 (90%)   = ', &
         real(n_in90,dp)/real(n_tot,dp), '   (1.0 => sin restriccion)'
    write(*,'(A,F8.4)')   '   frac Dchi2 < 2.300 (68%)   = ', &
         real(n_in68,dp)/real(n_tot,dp)
  end subroutine barrido_nsi

end program chi2_ideal_nest
