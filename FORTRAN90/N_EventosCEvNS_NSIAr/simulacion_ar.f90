!=======================================================================
! Programa: simulacion_ar (ROI Extendida 4-20 electrones)
!           Adaptado para Argón con 3 valores de EEE en paralelo.
!=======================================================================
program simulacion_ar
  use constants
  use quenching
  use xsections
  use flux, only: kopeikin_spectrum, mueller_spectrum, fission_frac, &
                  E_nu_min_flux, E_nu_max
  implicit none

  ! --- Configuración de la Región de Interés (ROI) para Argón ---
  integer, parameter :: bin_min = 4
  integer, parameter :: bin_max = 20

  integer, parameter :: n_T = 500, n_E = 2000
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, integrando
  real(dp) :: tasa_Comb, E_ee
  integer :: i_T, i_E, n_bin, i, j, k

  real(dp), parameter :: P_th = 3.1_dp, E_f_MeV = 204.0_dp, L_m = 19.0_dp
  real(dp), parameter :: dias_exposicion = 331.0_dp
  real(dp) :: flux_factor, atoms_per_kg, sec_per_day, conv_keV

  integer, parameter :: n_u = 1000, n_d = 1000
  real(dp) :: eps_y, eps_x, alpha_best
  real(dp) :: eps_y_min, eps_y_max, eps_x_min, eps_x_max, deps_y, deps_x
  real(dp), parameter :: sin2th_SM = 0.23857_dp, sigma_alpha = 0.169_dp
  real(dp) :: QW_SM, q_nsi_ee, q_nsi_emu, q_nsi_etau, q_eff2
  
  integer :: ipar
  integer :: u_out1, u_out2, u_out3, u_conf
  character(len=200) :: outdir, file_pes, file_nom, file_opt, configfile
  character(len=50)  :: xlabel, ylabel

  ! ================== ARREGLOS PARA LAS 3 EEE (Argón) ==================
  real(dp) :: val_eee(3) = [0.40_dp, 0.70_dp, 1.00_dp]
  
  ! --- Arrays con ROI extendida ---
  real(dp) :: R_unit_3(3, bin_min:bin_max), R_SM_3(3, bin_min:bin_max)
  real(dp) :: R_exp_3(3, bin_min:bin_max), sigma_exp_3(3, bin_min:bin_max)
  real(dp) :: R_th_3(3, bin_min:bin_max), chi2_3(3), sum1_3(3), sum2_3(3)

  ! ================== SELECCIÓN DEL CASO (1-15) ==================
  ipar = 5

  eps_y_min = -1.0_dp; eps_y_max =  1.0_dp
  eps_x_min = -1.0_dp; eps_x_max =  1.0_dp
  deps_y = (eps_y_max - eps_y_min) / (n_u - 1)
  deps_x = (eps_x_max - eps_x_min) / (n_d - 1)

  outdir = './' 
  file_pes = trim(outdir) // 'chi2_ar_pesimista.dat'
  file_nom = trim(outdir) // 'chi2_ar_nominal.dat'
  file_opt = trim(outdir) // 'chi2_ar_optimista.dat'
  configfile = trim(outdir) // 'nsi_configAr.txt'

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*sin2th_SM)/2.0_dp * Z_Ge
  QV2 = 1.0_dp 
  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day  = 86400.0_dp
  conv_keV     = 1000.0_dp
  flux_factor = (P_th * 1.0d9 / (E_f_MeV * 1.60218d-13)) / (4.0_dp * PI * (L_m * 100.0_dp)**2)

  ! ================== FASE 1: PRECALCULAR R_unit ==================
  R_unit_3 = 0.0_dp
  T_nr_min = 0.1_dp / 1000.0_dp; T_nr_max = 3.0_dp / 1000.0_dp
  dT = (T_nr_max - T_nr_min) / (n_T - 1); dT_keV = dT * 1000.0_dp
  dE = E_nu_max / (n_E - 1)

  write(*,*) 'Integrando espectro base para Argon (ROI 4-20)...'
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Comb = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle
        if (E_nu < E_nu_min_flux) then
           integrando = kopeikin_spectrum(E_nu) * dsigma_dT(E_nu, T_nr)
        else
           integrando = (fission_frac(1)*mueller_spectrum(E_nu,1) + &
                         fission_frac(2)*mueller_spectrum(E_nu,2) + &
                         fission_frac(3)*mueller_spectrum(E_nu,3) + &
                         fission_frac(4)*mueller_spectrum(E_nu,4)) * dsigma_dT(E_nu, T_nr)
        end if
        tasa_Comb = tasa_Comb + integrando * dE
     end do
     tasa_Comb = tasa_Comb * flux_factor * conv * atoms_per_kg * sec_per_day / conv_keV
     
     ! === MAPEO DINÁMICO DE UMBRALES (ROI extendida) ===
     E_ee = T_nr * QF(T_nr)
     do k = 1, 3
        n_bin = nint((E_ee / eta) * val_eee(k))
        if (n_bin >= bin_min .and. n_bin <= bin_max) then
           R_unit_3(k, n_bin) = R_unit_3(k, n_bin) + &
                (tasa_Comb * dT_keV * total_mass_kg * dias_exposicion)
        end if
     end do
  end do

  do k = 1, 3
     do n_bin = bin_min, bin_max
        R_SM_3(k, n_bin) = R_unit_3(k, n_bin) * (QW_SM**2)
        R_exp_3(k, n_bin) = R_SM_3(k, n_bin)
        sigma_exp_3(k, n_bin) = sqrt(R_SM_3(k, n_bin))
        if (sigma_exp_3(k, n_bin) < 1.0d-6) sigma_exp_3(k, n_bin) = 1.0d-6 
     end do
     write(*,*) 'EEE =', val_eee(k)
     write(*,*) '  Eventos en ROI:', sum(R_SM_3(k, bin_min:bin_max))
  end do

  ! ================== FASE 2: BARRIDO Y CHI CUADRADO ==================
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
  
  open(newunit=u_out1, file=file_pes, status='replace')
  open(newunit=u_out2, file=file_nom, status='replace')
  open(newunit=u_out3, file=file_opt, status='replace')
  
  write(*,*) 'Generando malla NSI para las 3 eficiencias de Argon (ROI 4-20)...'

  do i = 0, n_u-1
     eps_y = eps_y_min + i * deps_y 
     do j = 0, n_d-1
        eps_x = eps_x_min + j * deps_x 
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
        
        do k = 1, 3
           R_th_3(k, :) = R_unit_3(k, :) * q_eff2
           
           sum1_3(k) = sum((R_exp_3(k, :) * R_th_3(k, :)) / sigma_exp_3(k, :)**2)
           sum2_3(k) = sum(R_th_3(k, :)**2 / sigma_exp_3(k, :)**2)
           
           alpha_best = (sum1_3(k) - sum2_3(k)) / (sum2_3(k) + 1.0_dp/sigma_alpha**2)
           
           chi2_3(k) = sum(((R_exp_3(k, :) - (1.0_dp + alpha_best)*R_th_3(k, :))**2) / &
                       sigma_exp_3(k, :)**2) + (alpha_best/sigma_alpha)**2
        end do
        
        write(u_out1, '(3ES15.6)') eps_x, eps_y, chi2_3(1)
        write(u_out2, '(3ES15.6)') eps_x, eps_y, chi2_3(2)
        write(u_out3, '(3ES15.6)') eps_x, eps_y, chi2_3(3)
     end do
     write(u_out1, *)
     write(u_out2, *)
     write(u_out3, *)
  end do
  
  close(u_out1); close(u_out2); close(u_out3)
  write(*,*) '--- 3 Matrices NSI para Argon (ROI 4-20) Generadas con Exito ---'
end program simulacion_ar
