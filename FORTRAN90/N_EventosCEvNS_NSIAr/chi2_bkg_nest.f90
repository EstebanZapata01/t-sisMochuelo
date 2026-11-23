!=======================================================================
! Programa: chi2_bkg_nest (Ar)
!   Sensibilidad esperada con el fondo que declara RED-100, S/sqrt(B) ~ 4 a
!   62 kg*dia, repartido igual en N_e = 1..5. El fondo se supone conocido:
!       Dchi2(A) = (1-A)^2 W,   W = sum_k S_k^2/(S_k+B_k),
!       A_90 = 1 + sqrt(2.706/W)
!   Escenarios: (0) sin fondo; (1) con el fondo declarado.
!
!   Salidas (datos/): sensib_bkg_Ar.dat, bkg_bins_Ar.dat
!=======================================================================
program chi2_bkg_nest
  use constants
  use mod_tnr_to_e
  use xsections_nest
  use mod_stats,   only: binomial_prob
  use flux,        only: flujo_diferencial, E_nu_max
  implicit none

  integer, parameter :: NE_LO = 1, NE_HI = 5                 ! ROI
  integer,  parameter :: n_T = 800, n_E = 2000
  real(dp), parameter :: dchi2_90 = 2.706_dp
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, peso
  real(dp) :: tasa_Comb, QW_SM, p_F, dummy
  integer  :: i_T, i_E, k, n_F, u
  real(dp) :: S_bin(NE_LO:NE_HI)                             ! senal CEvNS SM [ev/(kg dia)]
  real(dp) :: expo_val, S_tot, B_tot, B_bin
  character(len=300) :: datadir

  datadir = '../../datos/'
  call inicializar_nest()
  dummy = flujo_diferencial(1.0_dp)

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2   = QW_SM**2

  ! ==================================================================
  ! 1. SENAL CEvNS SM  S(N_e)   [ev/(kg dia)]   (misma cadena que chi2_ideal)
  ! ==================================================================
  S_bin    = 0.0_dp
  T_nr_min = 0.20_dp / 1000.0_dp
  T_nr_max = 3.5_dp  / 1000.0_dp        ! T_max = 3.44 keV a E_nu = 8 MeV
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
  ! 2. FONDO: B_tot tal que S/sqrt(B) = SB_ref46, igual en cada bin
  ! ==================================================================
  expo_val = mass_Ar_kg * 1.0_dp
  S_tot    = expo_val * sum(S_bin)
  B_tot    = (S_tot / SB_ref46)**2
  B_bin    = B_tot / (expo_val * real(NE_HI - NE_LO + 1, dp))    ! ev/(kg dia) por bin

  write(*,'(/,A)')       ' ===== ARGON: sensibilidad esperada a 62 kg*dia ====='
  write(*,'(A,ES12.4)')  '   S_tot (eventos)                    = ', S_tot
  write(*,'(A,ES12.4)')  '   B_tot (S/sqrt(B) = 4, ref. [46])    = ', B_tot
  write(*,'(A,ES12.4,A)')'   B_k plano por bin                  = ', B_bin, ' ev/(kg dia)'
  write(*,'(A,F10.4)')   '   A_90 sin fondo                     = ', A90_fondo(0.0_dp)
  write(*,'(A,F10.4)')   '   A_90 con el fondo de la ref. [46]  = ', A90_fondo(B_bin)

  open(newunit=u, file=trim(datadir)//'sensib_bkg_Ar.dat', status='replace')
  write(u,'(A)') '# escenario   exposicion_kgd   A_90(xSM,Dchi2=2.706)   S_tot   B_tot   S/sqrtB'
  write(u,'(A12,ES16.6,F14.5,3ES14.5)') 'sin_fondo', expo_val, A90_fondo(0.0_dp), S_tot, 0.0_dp, 0.0_dp
  write(u,'(A12,ES16.6,F14.5,3ES14.5)') 'fondo_ref46', expo_val, A90_fondo(B_bin), S_tot, B_tot, SB_ref46
  close(u)

  open(newunit=u, file=trim(datadir)//'bkg_bins_Ar.dat', status='replace')
  write(u,'(A)') '# N_e  S_k  B_k   [ev/(kg dia)]; exposicion de la comparacion: 62 kg dia'
  do k = NE_LO, NE_HI
     write(u,'(I4,2ES16.7)') k, S_bin(k), B_bin
  end do
  close(u)

contains

  function A90_fondo(Bk) result(A90o)          ! Bk [ev/(kg dia)] por bin, igual en todos
    real(dp), intent(in) :: Bk
    real(dp) :: A90o, W, Sk
    integer  :: kk
    W = 0.0_dp
    do kk = NE_LO, NE_HI
       Sk = S_bin(kk) * expo_val
       if (Sk > 0.0_dp) W = W + Sk*Sk / (Sk + Bk * expo_val)
    end do
    A90o = 1.0_dp + sqrt(dchi2_90 / max(W, 1.0e-300_dp))
  end function A90_fondo

end program chi2_bkg_nest
