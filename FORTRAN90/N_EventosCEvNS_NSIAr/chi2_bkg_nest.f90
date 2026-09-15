!=======================================================================
! Programa: chi2_bkg_nest  (ARGON - sensibilidad esperada estilo RED-100 SV,
!                           con FONDO de 39Ar SIMULADO)
!
!   RED-100 SV: senal CEvNS SIMULADA + fondo -> Asimov -> A_90 en
!   Delta chi2 = 2.706.  RED-100 usa el fondo OFF MEDIDO.  Para argon nadie
!   ha medido nada, asi que el fondo dominante y especifico de argon -el
!   decaimiento beta del 39Ar (Q_beta = 565 keV, actividad publicada)- se
!   SIMULA sin inventar: isotopo conocido + forma de beta permitido de libro.
!
!   Todo en N_e (SIN paso PE: el PE es solo el formato del residuo publicado
!   de Xe; ver README_Ar.txt).
!
!   Senal   S(N_e): flujo x seccion eficaz x binomial(N_e), ROI N_e = 1..5.
!                   yield NR = tabla nest_Ar_218V_dense.txt (LArNEST anclado
!                   a ReD 2025; T_nr<2 keV = extrapolacion de modelo).
!   Fondo   B(N_e): espectro beta de 39Ar  dN/dT ~ F(Z,T) p_e E_e (Q-T)^2,
!                   plegado por el yield ER (tabla nest_Ar_ER_218V.txt) +
!                   binomial + EEE.  Normalizado por actividad x masa x tiempo.
!   Estadistica:    Delta chi2(A) = sum_k (1-A)^2 S_k^2/(S_k + B_k + s_SE*E)
!                   A_90 = 1 + sqrt(2.706 / sum_k S_k^2/(S_k+B_k+...))
!
!   EXPOSICION: UNA SOLA, real y justificada -- 62 kg (masa del "proximo
!   montaje" de ref.[46]) x 1 dia = 62 kg*dia, el MISMO numero que ref.[46]
!   declara para su propio S/sqrt(B)~4.  NO se toma prestada la exposicion
!   de Xe (192 kg*dia, sin ningun plan publicado de correr Ar esa cantidad
!   de tiempo) y NO se barre la exposicion (no hay ningun otro numero
!   publicado para barrer). El parametro que SI se mueve es el escenario de
!   FONDO -- un error/incertidumbre discreta y justificada, no la exposicion.
!
!   Cuatro escenarios discretos, todos a 62 kg*dia ("corridas de reemplazo",
!   Sec. brecha ideal->real de metodologia.tex, no una curva continua):
!     (0) sin fondo       -- techo intrinseco, reproduce chi2red100_nest.
!     (1) 39Ar UAr         -- actividad publicada (DarkSide-50, 7.3e-4 Bq/kg).
!     (2) 39Ar atmosferico -- actividad publicada (ref.[46], 1 Bq/kg).
!     (3) UAr + s_SE       -- (1) mas el piso de apilamiento de electron
!         unico CALIBRADO para reproducir el S/sqrt(B)~4 que ref.[46]
!         declara (ref.[46] deja ese fondo sin resolver; se infiere su
!         tamano implicito a partir del propio numero publicado, no se
!         supone nada mas). Este escenario SI se ejecuta (antes solo se
!         imprimia el s_SE necesario como diagnostico, sin volver a correr
!         el ajuste con el).
!
!   Observable fisico que se reporta (el mismo que usa ref.[46], no uno
!   inventado): S/sqrt(B) por escenario, mas A_90. Ambos se calculan SOLO a
!   partir de la prediccion SM (S_k) y del fondo simulado de actividades
!   publicadas (B_k) -- ninguna comparacion con datos de Xe.
!
!   Entrada : datos/nest_Ar_218V_dense.txt      (via mod_tnr_to_e, yield NR)
!             datos/nest_Ar_ER_218V.txt         (yield ER, lector propio)
!   Salida  : datos/sensib_bkg_Ar.dat
!             (escenario  exposicion_kgd(=62 fijo)  A_90  S_tot  B_tot  S/sqrtB)
!
! Tesis   : metodologia.tex Sec. "Resultado y validacion contra la ref. [46]".
!=======================================================================
program chi2_bkg_nest
  use constants
  use mod_tnr_to_e
  use xsections_nest
  use mod_stats,   only: binomial_prob
  use flux,        only: flujo_diferencial, E_nu_max
  implicit none

  ! ----- ROI -----
  integer, parameter :: NE_LO = 1, NE_HI = 5

  ! ----- grilla senal -----
  integer,  parameter :: n_T = 800, n_E = 2000
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, peso
  real(dp) :: tasa_Comb, QW_SM, p_F
  integer  :: i_T, i_E, k, n_F

  real(dp) :: S_bin(NE_LO:NE_HI)     ! tasa CEvNS SM  [ev/(kg dia)]

  ! ----- fondo 39Ar -----
  integer,  parameter :: n_ER = 6000
  real(dp), allocatable :: Eer_t(:), Qyer_t(:), Fer_t(:)
  integer  :: n_er_pts
  real(dp) :: B_shape(NE_LO:NE_HI)   ! fraccion de decaimientos con N_e_ext=k
  real(dp) :: B_rate_atm(NE_LO:NE_HI), B_rate_uar(NE_LO:NE_HI)  ! decaim/(kg dia)
  real(dp) :: Te, dTe, spec, norm_beta, lambda_er, p_Fer, Fer
  integer  :: i_er, n_Fer

  ! ----- estadistica -----
  real(dp), parameter :: me_keV = 510.99895_dp
  real(dp), parameter :: alpha_fs = 1.0_dp / 137.035999_dp
  real(dp), parameter :: dchi2_90 = 2.706_dp

  real(dp) :: A90, S_tot, B_tot, expo_val, sef_needed
  integer  :: im, u

  character(len=300) :: datadir, f_er, f_out
  character(len=32)  :: lbl
  real(dp) :: dummy

  datadir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  f_er    = trim(datadir)//'nest_Ar_ER_218V.txt'
  f_out   = trim(datadir)//'sensib_bkg_Ar.dat'

  call inicializar_nest()
  dummy = flujo_diferencial(1.0_dp)

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2   = QW_SM**2

  ! ==================================================================
  ! 1. SENAL CEvNS SM  S(N_e)   [ev/(kg dia)]   (misma cadena que chi2_ideal)
  ! ==================================================================
  S_bin    = 0.0_dp
  T_nr_min = 0.20_dp / 1000.0_dp
  T_nr_max = 3.5_dp  / 1000.0_dp        ! ref.[46] + cinematica E_nu=8 MeV
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
     tasa_Comb = tasa_Comb * (NA / A_Ge) * 1000.0_dp * 86400.0_dp
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)
     call obtener_nest_binomial(T_nr * 1000.0_dp, n_F, p_F)
     do k = NE_LO, NE_HI
        S_bin(k) = S_bin(k) + tasa_Comb * dT_keV * peso * &
                   binomial_prob(k, n_F, p_F * EEE)
     end do
  end do

  ! ==================================================================
  ! 2. FONDO 39Ar  B(N_e)   -- espectro beta permitido plegado por yield ER
  ! ==================================================================
  ! 2a. leer tabla ER (E_er[keV]  Qy_er[e-/keV]  F_er)
  allocate(Eer_t(n_ER), Qyer_t(n_ER), Fer_t(n_ER))
  call leer_tabla_ER(f_er, n_ER, Eer_t, Qyer_t, Fer_t, n_er_pts)

  ! 2b. integrar el espectro beta de 39Ar  (T de ~0 a Q_beta)
  B_shape   = 0.0_dp
  norm_beta = 0.0_dp
  dTe = Qbeta_Ar39_keV / real(n_ER - 1, dp)
  do i_er = 1, n_ER
     Te = (i_er - 1) * dTe
     if (Te <= 0.0_dp .or. Te >= Qbeta_Ar39_keV) cycle
     spec = beta_spectrum(Te, Qbeta_Ar39_keV, Z_daughter_Ar39)
     peso = merge(0.5_dp, 1.0_dp, i_er == 1 .or. i_er == n_ER)
     norm_beta = norm_beta + spec * peso * dTe

     ! <N_e> a esta energia de retroceso electronico
     lambda_er = Te * interp_lin(Eer_t, Qyer_t, n_er_pts, Te)
     if (lambda_er <= 0.0_dp) cycle
     Fer   = interp_lin(Eer_t, Fer_t, n_er_pts, Te)
     if (Fer < 0.02_dp) Fer = 0.02_dp
     if (Fer > 0.98_dp) Fer = 0.98_dp
     p_Fer = 1.0_dp - Fer
     n_Fer = nint(lambda_er / p_Fer)
     if (n_Fer < 1) n_Fer = 1
     do k = NE_LO, NE_HI
        B_shape(k) = B_shape(k) + spec * peso * dTe * &
                     binomial_prob(k, n_Fer, p_Fer * EEE)
     end do
  end do
  do k = NE_LO, NE_HI
     B_shape(k) = B_shape(k) / max(norm_beta, 1.0e-300_dp)   ! fraccion de decaim.
     ! decaimientos por kg y por dia con N_e_ext = k
     B_rate_atm(k) = act_Ar39_atm_Bq_kg * 86400.0_dp * B_shape(k)
     B_rate_uar(k) = act_Ar39_UAr_Bq_kg * 86400.0_dp * B_shape(k)
  end do

  ! ==================================================================
  ! 3. EXPOSICION UNICA Y JUSTIFICADA: 62 kg (ref.[46], "proximo montaje")
  !    x 1 dia = 62 kg*dia -- el mismo numero que ref.[46] declara para su
  !    S/sqrt(B)~4. NO se toma prestada la exposicion de Xe (192 kg*dia) y
  !    NO se barre la exposicion (no hay ningun otro numero publicado).
  ! ==================================================================
  expo_val = mass_Ar_kg * 1.0_dp

  write(*,'(/,A)') '===== SENSIBILIDAD ESPERADA (RED-100 SV) con FONDO 39Ar simulado - ARGON ====='
  write(*,'(A,I0,A,I0)')  ' ROI en N_e                            = ', NE_LO, ' .. ', NE_HI
  write(*,'(A,F8.4)')     ' EEE                                   = ', EEE
  write(*,'(A,F8.1,A)')   ' Exposicion (ref.[46]: 62 kg x 1 dia)  = ', expo_val, ' kg*dia'
  write(*,'(A)')          ' --- senal CEvNS SM por bin (ev/(kg dia)) ---'
  S_tot = 0.0_dp
  do k = NE_LO, NE_HI
     write(*,'(A,I2,A,ES13.5)') '   N_e = ', k, ' : ', S_bin(k)
     S_tot = S_tot + S_bin(k)
  end do
  write(*,'(A,ES13.5)')  '   S_tot                    = ', S_tot
  write(*,'(A)')          ' --- fondo 39Ar por bin (decaim/(kg dia)) [UAr | atmosferico] ---'
  do k = NE_LO, NE_HI
     write(*,'(A,I2,A,2ES13.5)') '   N_e = ', k, ' : ', B_rate_uar(k), B_rate_atm(k)
  end do

  ! ==================================================================
  ! 4. CALIBRAR el piso de apilamiento de electron unico (s_SE) contra el
  !    S/sqrt(B)~4 que ref.[46] declara a 62 kg*dia (UAr): es la unica
  !    forma de fijar ese parametro con un numero publicado, sin suponer
  !    nada mas. sef_needed = tamano IMPLICITO del fondo de apilamiento
  !    que ref.[46] deja sin resolver.
  ! ==================================================================
  block
    real(dp) :: B39_62, B_target
    integer  :: nb
    nb = NE_HI - NE_LO + 1
    call a90_con_fondo(1, expo_val, S_bin, B_rate_uar, B_rate_atm, &
                       NE_LO, NE_HI, 0.0_dp, A90, S_tot, B39_62)
    B_target   = (S_tot / SB_ref46)**2               ! B que da S/sqrt(B) = 4
    sef_needed = max((B_target - B39_62) / (expo_val * real(nb, dp)), 0.0_dp)
    write(*,'(/,A)') ' ===== CALIBRACION del piso s_SE contra ref.[46] (Physics 5, 492 (2023)) ====='
    write(*,'(A,ES12.4)')  '   B_tot (39Ar UAr solo, sin s_SE)   = ', B39_62
    write(*,'(A,F12.2)')   '   S / sqrt(B)  (solo 39Ar UAr)      = ', S_tot/sqrt(max(B39_62,1.0e-300_dp))
    write(*,'(A,F12.2)')   '   ref.[46] declara ~                = ', SB_ref46
    write(*,'(A)')         '   -> el 39Ar por solape espectral puro es MUCHO menor de lo que'
    write(*,'(A)')         '      sugiere ref.[46]: llegar a N_e<=5 pide E_er <~ 0.1 keV (cola'
    write(*,'(A)')         '      extrema del espectro beta + extrapolacion del yield ER).'
    write(*,'(A,ES12.4,A)')'   s_SE que reconcilia con ref.[46]  = ', sef_needed, ' ev/(kg dia) por bin'
    write(*,'(A)')         '      (= tamano implicito del apilamiento de SE que ref.[46] deja'
    write(*,'(A)')         '       SIN RESOLVER; se ejecuta como escenario (3) abajo.)'
  end block

  ! ==================================================================
  ! 5. CUATRO ESCENARIOS DISCRETOS, TODOS a la MISMA exposicion (62
  !    kg*dia): el "parametro de error" que se mueve es el fondo, no la
  !    exposicion ("corridas de reemplazo", no una curva continua).
  ! ==================================================================
  open(newunit=u, file=f_out, status='replace')
  write(u,'(A)') '# escenario                     exposicion_kgd   A_90(xSM,Dchi2=2.706)   S_tot   B_tot   S/sqrtB'
  write(*,'(/,A)') ' --- escenario                       A_90        S/sqrt(B) ---'

  do im = 0, 3
     select case(im)
     case(0)
        lbl = 'sin_fondo'
        call a90_con_fondo(0, expo_val, S_bin, B_rate_uar, B_rate_atm, &
                           NE_LO, NE_HI, 0.0_dp, A90, S_tot, B_tot)
     case(1)
        lbl = 'Ar39_UAr'
        call a90_con_fondo(1, expo_val, S_bin, B_rate_uar, B_rate_atm, &
                           NE_LO, NE_HI, 0.0_dp, A90, S_tot, B_tot)
     case(2)
        lbl = 'Ar39_atmosferico'
        call a90_con_fondo(2, expo_val, S_bin, B_rate_uar, B_rate_atm, &
                           NE_LO, NE_HI, 0.0_dp, A90, S_tot, B_tot)
     case(3)
        lbl = 'UAr_SEfloor_ref46'
        call a90_con_fondo(1, expo_val, S_bin, B_rate_uar, B_rate_atm, &
                           NE_LO, NE_HI, sef_needed, A90, S_tot, B_tot)
     end select
     write(*,'(A,A20,F12.4,ES14.4)') '   ', lbl, A90, &
          merge(S_tot/sqrt(max(B_tot,1.0e-300_dp)), -1.0_dp, im /= 0)
     write(u,'(A24,ES16.6,F14.5,3ES14.5)') adjustl(lbl), expo_val, A90, &
          S_tot, B_tot, merge(S_tot/sqrt(max(B_tot,1.0e-300_dp)), 0.0_dp, im /= 0)
  end do
  close(u)
  write(*,'(/,A)') '  Salida: datos/sensib_bkg_Ar.dat (4 escenarios, 62 kg*dia fijo)'

  deallocate(Eer_t, Qyer_t, Fer_t)

contains

  !---------------------------------------------------------------------
  ! espectro beta permitido de 39Ar:  dN/dT ~ F(Z,T) p_e E_e (Q - T)^2
  ! T = energia cinetica del electron [keV]; Fermi no-relativista.
  !---------------------------------------------------------------------
  function beta_spectrum(T, Q, Zd) result(s)
    real(dp), intent(in) :: T, Q, Zd
    real(dp) :: s, Ee, pe, eta, Ffermi
    if (T <= 0.0_dp .or. T >= Q) then
       s = 0.0_dp; return
    end if
    Ee = T + me_keV
    pe = sqrt(T*(T + 2.0_dp*me_keV))
    eta = alpha_fs * Zd * Ee / pe            ! beta-: atraccion (+)
    Ffermi = 2.0_dp*pi*eta / (1.0_dp - exp(-2.0_dp*pi*eta))
    s = Ffermi * pe * Ee * (Q - T)**2
  end function beta_spectrum

  !---------------------------------------------------------------------
  function interp_lin(x, y, n, xq) result(yq)
    real(dp), intent(in) :: x(:), y(:), xq
    integer,  intent(in) :: n
    real(dp) :: yq, f
    integer  :: lo, hi, mid
    if (xq <= x(1)) then
       yq = y(1); return
    else if (xq >= x(n)) then
       yq = y(n); return
    end if
    lo = 1; hi = n
    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       if (xq >= x(mid)) then; lo = mid; else; hi = mid; end if
    end do
    f = (xq - x(lo)) / (x(hi) - x(lo))
    yq = y(lo) + f * (y(hi) - y(lo))
  end function interp_lin

  !---------------------------------------------------------------------
  subroutine leer_tabla_ER(fname, nmax, xe, qy, fe, npts)
    character(len=*), intent(in)  :: fname
    integer,          intent(in)  :: nmax
    real(dp),         intent(out) :: xe(nmax), qy(nmax), fe(nmax)
    integer,          intent(out) :: npts
    integer :: uu, ios
    character(len=300) :: linea
    real(dp) :: a, b, c
    open(newunit=uu, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: no se encontro ', trim(fname)
       write(*,*) '  Genera datos/nest_Ar_ER_218V.txt con  python/nest_Ar.py'
       stop 1
    end if
    npts = 0
    do
       read(uu, '(A)', iostat=ios) linea
       if (ios /= 0) exit
       linea = adjustl(linea)
       if (linea(1:1) == '#' .or. len_trim(linea) == 0) cycle
       read(linea, *, iostat=ios) a, b, c
       if (ios /= 0) cycle
       npts = npts + 1
       if (npts > nmax) then
          write(*,*) 'ERROR: tabla ER mas larga que nmax=', nmax
          stop 1
       end if
       xe(npts) = a; qy(npts) = b; fe(npts) = c
    end do
    close(uu)
    write(*,'(A,I0,A)') ' Tabla ER 39Ar: ', npts, ' puntos leidos.'
  end subroutine leer_tabla_ER

  !---------------------------------------------------------------------
  ! A_90 con fondo:  W = sum_k S_k^2 / (S_k + B_k + se_floor*E)
  !   escenario im: 0 = sin fondo, 1 = UAr, 2 = atmosferico
  !---------------------------------------------------------------------
  subroutine a90_con_fondo(im, E, Sb, Bu, Ba, klo, khi, sef, A90o, Sto, Bto)
    integer,  intent(in)  :: im, klo, khi
    real(dp), intent(in)  :: E, Sb(klo:khi), Bu(klo:khi), Ba(klo:khi), sef
    real(dp), intent(out) :: A90o, Sto, Bto
    real(dp) :: Sk, Bk, Wk
    integer  :: kk
    Wk = 0.0_dp; Sto = 0.0_dp; Bto = 0.0_dp
    do kk = klo, khi
       Sk = Sb(kk) * E
       select case (im)
       case (0); Bk = 0.0_dp
       case (1); Bk = Bu(kk) * E
       case (2); Bk = Ba(kk) * E
       end select
       Bk = Bk + sef * E
       Sto = Sto + Sk
       Bto = Bto + Bk
       if (Sk > 0.0_dp) Wk = Wk + Sk*Sk / (Sk + Bk)
    end do
    A90o = 1.0_dp + sqrt(dchi2_90 / max(Wk, 1.0e-300_dp))
  end subroutine a90_con_fondo

end program chi2_bkg_nest
