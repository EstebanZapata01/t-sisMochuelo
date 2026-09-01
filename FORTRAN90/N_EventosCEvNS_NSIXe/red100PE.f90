!=======================================================================
! Programa: red100PE_detallado (Plantillas PCHIP - Flujo Explícito)
! Rol     : plantilla del espectro CE$\nu$NS SM en energia corregida [PE]
!           (convoluciona N_e con SEG y sig1 de mod_detector).
! Pipeline: etapa "N_e -> PE" -> escribe
!           datos/ionization_spectra_detallado.dat, que es la entrada de
!           chi2.f90 (ajuste 1D ON-OFF).
! Tesis   : metodologia.tex Sec. 3.3 y Sec. 4.
! Decision metodologica clave: la validacion compara tasa_ion_extraidos
!   contra la Fig. 3 y contra las dos curvas digitalizadas de la Fig. 6
!   (retencion global sum(f6_a)/sum(f6_b) ~ 0.2555). En la carpeta de Ar
!   este espectro en PE NO es fisico (mod_detector es de LXe).
!=======================================================================
program red100PE_detallado
  use constants
  use mod_tnr_to_e        
  use xsections_nest      
  use mod_stats           
  use mod_detector        
  use flux, only: flujo_diferencial, E_nu_max
  implicit none

  integer, parameter :: n_T = 500, n_E = 2000, n_ion = 15
  real(dp) :: T_nr, T_nr_min, T_nr_max, dT, dT_keV, E_nu, dE, peso
  integer :: i_T, i_E, i_bin, u_out, i_pe
  real(dp) :: tasa_Comb, S_pe, bin_width_pe, QW_SM
  real(dp), allocatable :: array_tasa_Comb(:), tasa_ion_extraidos(:), contribuciones(:)
  real(dp) :: total_bin
  integer  :: n_F, ipk
  real(dp) :: p_F
  character(len=250) :: outdir, filename, datafile
  real(dp) :: atoms_per_kg, sec_per_day, dummy
  real(dp) :: R_tot_comb, sum_ext_all, sum_ext_roi, sum_pe_total, pe_pk_val

  bin_width_pe = 5.0_dp
  
  outdir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  datafile = trim(outdir) // 'templates_SE_pchip.dat'
  filename = trim(outdir) // 'ionization_spectra_detallado.dat'

  call inicializar_nest()
  call cargar_SE_data(datafile)
  
  ! Disparamos la inicialización automática del flujo
  dummy = flujo_diferencial(1.0_dp)

  QW_SM = -N_Ge/2.0_dp + (1.0_dp - 4.0_dp*s2w)/2.0_dp * Z_Ge
  QV2 = QW_SM**2

  atoms_per_kg = (NA / A_Ge) * 1000.0_dp
  sec_per_day = 86400.0_dp

  allocate(array_tasa_Comb(n_T), tasa_ion_extraidos(n_ion), contribuciones(7))
  tasa_ion_extraidos(:) = 0.0_dp

  T_nr_min = 0.20_dp / 1000.0_dp; T_nr_max = 3.0_dp / 1000.0_dp   ! 0.20 keV = suelo de NEST v2.4.0
  dT = (T_nr_max - T_nr_min) / (n_T - 1); dT_keV = dT * 1000.0_dp
  dE = E_nu_max / (n_E - 1)

  write(*,*) "-> Integrando el espectro continuo..."
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     tasa_Comb = 0.0_dp
     do i_E = 1, n_E
        E_nu = (i_E - 1) * dE
        if (E_nu < sqrt(M_Ge * T_nr / 2.0_dp)) cycle
        peso = merge(0.5_dp, 1.0_dp, i_E == 1 .or. i_E == n_E)
        tasa_Comb = tasa_Comb + (flujo_diferencial(E_nu) * dsigma_dT(E_nu, T_nr) * dE * peso)
     end do
     array_tasa_Comb(i_T) = tasa_Comb * atoms_per_kg * sec_per_day
  end do

  R_tot_comb = 0.0_dp
  do i_T = 1, n_T
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)
     R_tot_comb = R_tot_comb + array_tasa_Comb(i_T) * dT_keV * peso
  end do

  write(*,*) "-> Aplicando fluctuacion de N_e (Binomial sub-Poissoniana de NEST)..."
  do i_T = 1, n_T
     T_nr = T_nr_min + (i_T - 1) * dT
     call obtener_nest_binomial(T_nr * 1000.0_dp, n_F, p_F)
     peso = merge(0.5_dp, 1.0_dp, i_T == 1 .or. i_T == n_T)

     ! Espectro TEORICO (Fig. 3 del paper): sin cortes de seleccion (eff_ROI) y
     ! sin livetime_frac. N_e extraidos ~ Binomial(n_F, p_F*EEE).
     do i_bin = 1, n_ion
        tasa_ion_extraidos(i_bin) = tasa_ion_extraidos(i_bin) + &
     (array_tasa_Comb(i_T) * dT_keV * peso * binomial_prob(i_bin, n_F, p_F * EEE))
     end do
  end do


  open(newunit=u_out, file=filename, status='replace')
  write(u_out, '(A)') '# PE_center   Total   1SE   2SE   3SE   4SE   5SE   6SE   7SE'

  write(*,*) "-> Escribiendo formato PCHIP..."
  sum_pe_total = 0.0_dp; pe_pk_val = 0.0_dp; ipk = 1
  do i_pe = 1, 400
     S_pe = (real(i_pe, dp) - 0.5_dp) * bin_width_pe
     total_bin = 0.0_dp

     do i_bin = 1, 7
        contribuciones(i_bin) = tasa_ion_extraidos(i_bin) * respuesta_empirica(S_pe, i_bin) * bin_width_pe
        total_bin = total_bin + contribuciones(i_bin)
     end do
     sum_pe_total = sum_pe_total + total_bin
     if (total_bin > pe_pk_val) then; pe_pk_val = total_bin; ipk = i_pe; end if

     write(u_out, '(F8.1, E15.6)', advance='no') S_pe, total_bin
     do i_bin = 1, 7
        write(u_out, '(E15.6)', advance='no') contribuciones(i_bin)
     end do
     write(u_out, *)
  end do
  close(u_out)

  sum_ext_all = sum(tasa_ion_extraidos(1:n_ion))
  sum_ext_roi = tasa_ion_extraidos(4) + tasa_ion_extraidos(5) &
              + tasa_ion_extraidos(6) + tasa_ion_extraidos(7)

  write(*,'(/,A)') '============ VALIDACION red100PE (vs arXiv:2411.18641 Figs. 3 y 6) ============'
  write(*,'(A,ES13.5,A)') ' Tasa CEvNS integrada 0.1-3 keV          = ', R_tot_comb, ' eventos/(kg dia)'
  write(*,'(A,ES13.5)')   ' Sum tasa_ion_extraidos(1:15) [eventos]  = ', sum_ext_all
  write(*,'(A,F8.4,A)')   '   fraccion con >= 1 e- extraido         = ', sum_ext_all/max(R_tot_comb,1.0e-30_dp), &
       '   (<<1: casi todo el retroceso CEvNS esta bajo el umbral de NEST 0.2 keV)'
  write(*,'(/,A)') ' N_e extraidos por bin [eventos/(kg dia)]  (comparar Fig. 3 top, "after extraction"):'
  do i_bin = 1, 7
     write(*,'(A,I2,A,ES13.5)') '   Ne = ', i_bin, ' : ', tasa_ion_extraidos(i_bin)
  end do
  write(*,'(A,ES13.5)') '   Sum ROI (Ne 4-7, antes de eff_ROI)   = ', sum_ext_roi

  ! ---- Comparacion bin a bin contra las figuras digitalizadas del paper ----
  ! pap  = Fig. 3  "N_e extraidos"                    (sin cortes de seleccion)
  ! f6_b = Fig. 6 (abajo) "CEvNS signal, before cuts" ["wpd_datasets (4).csv"]
  ! f6_a = Fig. 6 (abajo) "CEvNS signal, after cuts"  ["wpd_datasets (4).csv"]
  ! Todo en eventos/(kg dia), Ne = 4,5,6,7. eff_ROI(k) = f6_a(k)/f6_b(k).
  block
    real(dp) :: pap(4), f6_b(4), f6_a(4), sim(4)
    integer  :: b
    pap  = (/ 0.030182192_dp, 0.004289422_dp, 0.000575605_dp, 0.0000772415_dp /)
    f6_b = (/ 0.014513251_dp, 0.012017993_dp, 0.002075998_dp, 0.000253140_dp /)
    f6_a = (/ 0.001987556_dp, 0.003931315_dp, 0.001267480_dp, 0.000186642_dp /)
    sim = (/ tasa_ion_extraidos(4), tasa_ion_extraidos(5), &
             tasa_ion_extraidos(6), tasa_ion_extraidos(7) /)
    write(*,'(/,A)') ' Comparacion bin a bin  (sim = tasa_ion_extraidos, SIN eff_ROI):'
    write(*,'(A)')   '   Ne     sim           Fig.3         Fig.6 before   sim/Fig.3'
    do b = 1, 4
       write(*,'(I5,3ES15.5,F11.3)') b+3, sim(b), pap(b), f6_b(b), sim(b)/pap(b)
    end do
    write(*,'(A)') '   -- razones entre bins consecutivos (test de forma) --'
    write(*,'(A)') '   Ne(k)->k+1     sim          Fig.3        Fig.6 before'
    do b = 1, 3
       write(*,'(I8,A,I1,3F13.3)') b+3, '->', b+4, sim(b)/sim(b+1), &
            pap(b)/pap(b+1), f6_b(b)/f6_b(b+1)
    end do
    write(*,'(A)') '   -- senal en ROI tras cortes:  sim*eff_ROI  vs  Fig.6 after --'
    write(*,'(A)') '   Ne     sim*eff_ROI    Fig.6 after    (sim*eff)/Fig.6a'
    do b = 1, 4
       write(*,'(I5,2ES15.5,F13.3)') b+3, sim(b)*eff_ROI(b+3), f6_a(b), &
            sim(b)*eff_ROI(b+3)/f6_a(b)
    end do
    write(*,'(A,F7.4)') '   Retencion global Fig. 6   sum(after)/sum(before)      = ', &
         sum(f6_a) / sum(f6_b)
    write(*,'(A,F7.4)') '   Cociente de normalizacion sum(sim*eff_ROI)/sum(Fig.6a) = ', &
         sum(sim*eff_ROI(4:7)) / sum(f6_a)
    write(*,'(A)')      '   (0.2555 = el "75% signal loss in ROI" del texto; es cociente'
    write(*,'(A)')      '    de las dos curvas de la Fig. 6, no un promedio de eff por bin.)'
  end block
  write(*,'(/,A)') ' Espectro en PE:'
  write(*,'(A,F7.1,A,ES13.5)') '   pico en PE = ', (real(ipk,dp)-0.5_dp)*bin_width_pe, ' , valor = ', pe_pk_val
  write(*,'(A,ES13.5)') '   Sum_PE(Total) sobre todo el espectro = ', sum_pe_total
  write(*,'(A,F8.4,A)') '   Sum_PE(Total) / Sum tasa_ion(1:7)    = ', &
       sum_pe_total / max(tasa_ion_extraidos(1)+tasa_ion_extraidos(2)+tasa_ion_extraidos(3) &
                        + sum_ext_roi, 1.0e-30_dp), '   (~1 si las gaussianas caben en 0-2000 PE)'
  write(*,'(A,/)') '======================================================================================='

  write(*,*) "=== EXITOSO: espectro PE teorico (1SE..7SE), Eventos/(5PE*kg*dia) ==="
end program red100PE_detallado
