!=======================================================================
! Programa: chi2red100_nest  (ARGÓN - matriz NSI 2D)
!
!   Clon de ../N_EventosCEvNS_NSIXe/chi2red100_nest.f90 adaptado a Ar-40.
!
!   DIFERENCIA CLAVE con la versión de Xe: para Ar NO hay datos ON-OFF
!   reales (ni fondo OFF medido ni fluctuación ON). Por eso la estadística
!   es una SENSIBILIDAD ASIMOV DE CONTEO PURO en el ROI:
!
!       dN_k = R_k            (Asimov: el dato es la señal SM esperada)
!       sigma_k = sqrt(N_k)   (Poisson de SEÑAL; SIN fondo, SIN sistemáticos)
!       chi2(A) = sum_k (N_k - A*N_k)^2 / N_k = (1 - A)^2 * N_tot_ROI
!       A_best = 1,  chi2_min = 0
!
!   La NSI entra solo por la carga efectiva:
!       q_eff2(eps) = (Q_W + q_ee)^2 + q_emu^2 + q_etau^2
!       A_amp(eps)  = q_eff2 / Q_W^2      (A=1 -> SM)
!       Delta chi2(eps) = (1 - A_amp)^2 * N_tot_ROI
!
!   IDEALIZACIÓN (metodología acordada): esta rama mide la RESPUESTA
!   INTRÍNSECA del blanco (sigma_CEvNS . N^2 . charge yield . umbral), NO la
!   sensibilidad experimental alcanzable. Es la única fórmula que se puede
!   aplicar de forma idéntica y honesta a Xe y a Ar (ninguno usa fondo real
!   aquí). El límite OBSERVADO real de Xe (~111 ×SM) se reporta aparte y NO
!   se compara directamente con este número.
!   Para Ar solo es físico si se asume ARGÓN DEPLETADO (UAr): con Ar
!   atmosférico el fondo de 39Ar (~1 Bq/kg) domina el ROI por ~1e6-1e7.
!
!   Entrada : nada (calcula el espectro SM internamente, como mainred100_nest)
!             lee datos/nest_Ar_218V_dense.txt vía mod_tnr_to_e
!   Salida  : datos/chi2_nsi_2DAr.dat   (eps_x  eps_y  Delta_chi2)
!             datos/nsi_configAr.txt    (etiquetas de los ejes)
!             datos/sensib_expo_Ar.dat  (multiplicador  exposicion_kgd  A_90)
!
!   Caso por defecto ipar=5: eps_ee^dV (x) vs eps_emu^dV (y)
!
! Pipeline: variante "ventana fisica de Ar" (N_e = 1..5, T_nr_max=6 keV).
! Tesis   : metodologia.tex Sec. 5 (formalismo NSI), Sec. 6 (CONUS+ vs
!           RED-100) y Sec. 9. El mapeo eps -> q_nsi (15 casos ipar) es
!           IDENTICO a 2pchi2.f90 de CONUS+.
! Decision metodologica clave: A_amp(eps) = q_eff2/Q_w^2; sin datos de Ar,
!   la estadistica es Asimov de conteo puro (no residuo ON-OFF).
!=======================================================================
program chi2red100_Ar
  use constants
  use mod_tnr_to_e
  use xsections_nest
  use mod_stats,   only: binomial_prob
  use flux,        only: flujo_diferencial, E_nu_max
  implicit none

  ! ----- grilla de integración -----
  integer,  parameter :: n_T = 500, n_E = 2000
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, peso
  real(dp) :: tasa_Comb, QW_SM
  integer  :: i_T, i_E, n_bin, n_F
  real(dp) :: p_F

  ! ----- ROI en N_e (idealización común: eff_ROI = 1) -----
  ! ROI de Ar según §VII de arXiv:2411.18641: "below five ionization electrons"
  ! (ref. [46]) -> N_e = 1..5. (En Xe era 4..7, definido por el exceso de SE;
  ! §VII dice que en Ar el fondo de SE es "unclear", por eso conviene correr
  ! también 2..5 como variante — cambiar NE_LO a 2.)
  integer,  parameter :: NE_LO = 1, NE_HI = 5
  real(dp) :: R_roi(NE_LO:NE_HI)     ! tasa SM por bin  [eventos/(kg dia)]
  real(dp) :: N_roi(NE_LO:NE_HI)     ! cuentas SM absolutas (= R_roi * exposición)
  real(dp) :: Ntot_roi

  ! ----- NSI 2D -----
  integer,  parameter :: ipar = 5           ! eps_ee^dV (x)  vs  eps_emu^dV (y)
  integer,  parameter :: n_u = 1000, n_d = 1000
  real(dp) :: eps_x, eps_y, eps_min, eps_max, deps
  real(dp) :: q_nsi_ee, q_nsi_emu, q_nsi_etau, q_eff2, A_amp, dchi2
  real(dp) :: Aeff_lo, Aeff_hi
  integer  :: i, j, u_out, u_conf, u_sens, n_in90, n_in68, n_tot

  ! ----- sensibilidad Asimov vs exposición -----
  real(dp) :: A90
  real(dp), parameter :: mult(7) = &
       (/ 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 50.0_dp, 100.0_dp, 335.0_dp /)
  ! 335 x 192 kg*dia ~ 1 año astronómico a ~126 kg (comparar 15-20 ×SM, §VII)

  character(len=250) :: datadir, f_out, f_conf, f_sens
  character(len=60)  :: xlabel, ylabel
  real(dp) :: dummy

  datadir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  f_out   = trim(datadir)//'chi2_nsi_2DAr.dat'
  f_conf  = trim(datadir)//'nsi_configAr.txt'
  f_sens  = trim(datadir)//'sensib_expo_Ar.dat'

  call inicializar_nest()
  dummy = flujo_diferencial(1.0_dp)         ! dispara init del flujo

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2   = QW_SM**2

  ! ==================================================================
  ! 1. ESPECTRO SM EN LA ROI (bins de N_e)  -- igual que mainred100_nest
  ! ==================================================================
  R_roi = 0.0_dp
  T_nr_min = 0.20_dp / 1000.0_dp
  ! Ar-40: T_max = 2 E_nu^2/(M + 2 E_nu) ~ 3.44 keV para E_nu = 8 MeV y
  ! ~5.37 keV para E_nu = 10 MeV (= E_nu_max del flujo). 6.0 keV cubre con
  ! margen; la tabla nest_Ar_218V_dense.txt tambien llega a 6 keV.
  T_nr_max = 3.5_dp / 1000.0_dp   ! ref.[46] + cinematica E_nu=8 MeV (3.44 keV); antes 6.0
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
     tasa_Comb = tasa_Comb * (NA / A_Ge) * 1000.0_dp * 86400.0_dp   ! ev/(kg dia keV)
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)

     call obtener_nest_binomial(T_nr * 1000.0_dp, n_F, p_F)
     do n_bin = NE_LO, NE_HI
        R_roi(n_bin) = R_roi(n_bin) + tasa_Comb * dT_keV * peso * &
                       eff_ROI(n_bin) * binomial_prob(n_bin, n_F, p_F * EEE)
     end do
  end do

  Ntot_roi = 0.0_dp
  do n_bin = NE_LO, NE_HI
     N_roi(n_bin) = R_roi(n_bin) * exposure_ON_kgd
     Ntot_roi = Ntot_roi + N_roi(n_bin)
  end do

  ! ==================================================================
  ! 2. SENSIBILIDAD ASIMOV DE CONTEO PURO
  !    chi2(A) = (1-A)^2 * N_tot_ROI  ;  A_90 (1 gdl) : Dchi2 = 2.706
  ! ==================================================================
  write(*,'(/,A)') '========== SENSIBILIDAD Asimov conteo puro (ARGÓN) =========='
  write(*,'(A,F12.5)')  ' Q_w (Ar: N=22, Z=18)          = ', QW_SM
  write(*,'(A,F12.3)')  ' Q_w^2                          = ', QW_SM**2
  write(*,'(A,F10.3,A)')' EEE (Ar, ReD/DarkSide-50)     = ', EEE, ''
  write(*,'(A,I0,A,I0)')' ROI en N_e                     = ', NE_LO, ' .. ', NE_HI
  write(*,'(A)')        ' Tasa SM por bin  [eventos/(kg dia)] :'
  do n_bin = NE_LO, NE_HI
     write(*,'(A,I2,A,ES13.5)') '   N_e=', n_bin, ' : ', R_roi(n_bin)
  end do
  write(*,'(A,ES13.5)') ' Tasa total ROI   [eventos/(kg dia)] = ', sum(R_roi)
  write(*,'(A,F10.1,A)')' Exposición base                = ', exposure_ON_kgd, ' kg*dia'
  write(*,'(A,ES13.5)') ' N_tot_ROI (Asimov, base)       = ', Ntot_roi
  if (Ntot_roi > 0.0_dp) then
     write(*,'(A,F10.3,A)') ' A_90 (base, Dchi2=2.706)      = ', &
          1.0_dp + sqrt(2.706_dp / Ntot_roi), ' x SM'
  end if

  ! ---- sensibilidad vs exposición ----
  open(newunit=u_sens, file=f_sens, status='replace')
  write(u_sens,'(A)') '# multiplicador   exposicion_kgd   A_90(xSM, Dchi2=2.706, Asimov conteo puro)'
  write(*,'(/,A)') ' A_90 vs exposición (Asimov conteo puro, sin fondo):'
  write(*,'(A)')   '   x_base     kg*dia        A_90[xSM]'
  do i = 1, size(mult)
     if (Ntot_roi > 0.0_dp) then
        A90 = 1.0_dp + sqrt(2.706_dp / (Ntot_roi * mult(i)))
     else
        A90 = -1.0_dp
     end if
     write(*,'(F8.1,ES14.5,F14.3)') mult(i), exposure_ON_kgd*mult(i), A90
     write(u_sens,'(F10.1,2ES16.6)') mult(i), exposure_ON_kgd*mult(i), A90
  end do
  close(u_sens)

  ! ==================================================================
  ! 3. ETIQUETAS DE EJES (solo ipar=5; para otros casos ver la versión Xe)
  ! ==================================================================
  select case (ipar)
  case (5); xlabel = '$\epsilon_{ee}^{dV}$'; ylabel = '$\epsilon_{e\mu}^{dV}$'
  case default; xlabel = 'eps_x'; ylabel = 'eps_y'
  end select
  open(newunit=u_conf, file=f_conf, status='replace')
  write(u_conf,'(A)') trim(xlabel); write(u_conf,'(A)') trim(ylabel)
  close(u_conf)

  ! ==================================================================
  ! 4. BARRIDO NSI 2D   (Delta chi2 = (1 - A_amp)^2 * N_tot_ROI, exposición base)
  ! ==================================================================
  eps_min = -1.0_dp; eps_max = 1.0_dp
  deps = (eps_max - eps_min) / real(n_u - 1, dp)

  open(newunit=u_out, file=f_out, status='replace')
  write(u_out,'(A)') '# eps_x   eps_y   Delta_chi2   (Asimov conteo puro, exposicion base)'

  n_in90 = 0; n_in68 = 0; n_tot = 0
  Aeff_lo = 1.0e30_dp; Aeff_hi = -1.0e30_dp

  do i = 0, n_u-1
     eps_y = eps_min + i * deps
     do j = 0, n_d-1
        eps_x = eps_min + j * deps
        q_nsi_ee = 0.0_dp; q_nsi_emu = 0.0_dp; q_nsi_etau = 0.0_dp

        ! ipar = 5 : eps_ee^dV (x)  y  eps_emu^dV (y)   (mismas fórmulas que Xe)
        q_nsi_ee  = (eps_x)*Z_Ge + (2.0_dp*eps_x)*N_Ge
        q_nsi_emu = (eps_y)*Z_Ge + (2.0_dp*eps_y)*N_Ge

        q_eff2 = (QW_SM + q_nsi_ee)**2 + q_nsi_emu**2 + q_nsi_etau**2
        A_amp  = q_eff2 / QW_SM**2

        dchi2 = (1.0_dp - A_amp)**2 * Ntot_roi

        n_tot = n_tot + 1
        if (dchi2 <= 4.605_dp) n_in90 = n_in90 + 1     ! 90% 2 gdl
        if (dchi2 <= 2.300_dp) n_in68 = n_in68 + 1     ! 68% 2 gdl
        if (A_amp < Aeff_lo) Aeff_lo = A_amp
        if (A_amp > Aeff_hi) Aeff_hi = A_amp

        write(u_out,'(3ES15.6)') eps_x, eps_y, dchi2
     end do
     write(u_out,*)
  end do
  close(u_out)

  write(*,'(/,A)') ' --- Barrido NSI Ar (ipar=5: eps_ee^dV vs eps_emu^dV, |eps| <= 1) ---'
  write(*,'(A,2F10.3)') ' A_amp = q_eff2/Q_w^2 en la grilla: min, max = ', Aeff_lo, Aeff_hi
  write(*,'(A,F8.4,A)') ' Fracción del plano con Dchi2 < 4.605 (90%, 2 gdl) = ', &
       real(n_in90,dp)/real(n_tot,dp), '   (1.0 => sin restricción)'
  write(*,'(A,F8.4)')   ' Fracción del plano con Dchi2 < 2.300 (68%, 2 gdl) = ', &
       real(n_in68,dp)/real(n_tot,dp)
  write(*,'(A)') '================================================================================'

end program chi2red100_Ar
