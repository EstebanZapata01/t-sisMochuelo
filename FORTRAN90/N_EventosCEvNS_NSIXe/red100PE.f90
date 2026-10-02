! Plantilla del espectro CEvNS SM en energia corregida [PE] (convoluciona
! N_e con SEG y sig1 de mod_detector) -> datos/ionization_spectra_detallado.dat,
! la entrada de chi2.f90. Dos espacios que no se mezclan: N_e VERDADERO
! (tasa_ion_extraidos) y N_e RECONSTRUIDO (PE/27, ventanas +-0,5 e-
! recortadas a 110-189 PE, via prob_migracion); el total es la suma de las
! plantillas k=1..15, incluida la migracion desde k<4. Este espectro en PE
! no es fisico en la carpeta de Ar (mod_detector es especifico de LXe).
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
  real(dp) :: total_bin, c_k
  integer  :: n_F, ipk
  real(dp) :: p_F
  character(len=250) :: outdir, filename, datafile
  real(dp) :: atoms_per_kg, sec_per_day, dummy
  real(dp) :: R_tot_comb, sum_ext_all, sum_ext_roi, sum_pe_total, pe_pk_val

  bin_width_pe = 5.0_dp
  
  outdir = '../../datos/'
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
  write(u_out, '(A)') '# PE_center   Total(1SE..15SE)   1SE   2SE   3SE   4SE   5SE   6SE   7SE'

  write(*,*) "-> Escribiendo formato PCHIP..."
  sum_pe_total = 0.0_dp; pe_pk_val = 0.0_dp; ipk = 1
  do i_pe = 1, 400
     S_pe = (real(i_pe, dp) - 0.5_dp) * bin_width_pe
     total_bin = 0.0_dp

     do i_bin = 1, n_ion
        c_k = tasa_ion_extraidos(i_bin) * respuesta_gauss(S_pe, i_bin) * bin_width_pe
        total_bin = total_bin + c_k
        if (i_bin <= 7) contribuciones(i_bin) = c_k
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

  ! ---- Comparacion contra las figuras digitalizadas del paper ----
  ! pap  = Fig. 3  "N_e extraidos"                    (N_e VERDADERO)
  ! f6_b = Fig. 6 (abajo) "CEvNS signal, before cuts" (N_e RECONSTRUIDO)
  ! f6_a = Fig. 6 (abajo) "CEvNS signal, after cuts"  (N_e RECONSTRUIDO)
  ! Todo en eventos/(kg dia), Ne = 4,5,6,7. eff_ROI(j) = f6_a(j)/f6_b(j).
  block
    real(dp) :: pap(4), f6_b(4), f6_a(4), sim(4), rec(4)
    integer  :: b, k_e
    pap  = (/ 0.030182192_dp, 0.004289422_dp, 0.000575605_dp, 0.0000772415_dp /)
    f6_b = (/ 0.014513251_dp, 0.012017993_dp, 0.002075998_dp, 0.000253140_dp /)
    f6_a = (/ 0.001987556_dp, 0.003931315_dp, 0.001267480_dp, 0.000186642_dp /)
    sim = (/ tasa_ion_extraidos(4), tasa_ion_extraidos(5), &
             tasa_ion_extraidos(6), tasa_ion_extraidos(7) /)
    do b = 1, 4     ! N_e reconstruido: sum_k R_k * P(k -> j)
       rec(b) = 0.0_dp
       do k_e = 1, n_ion
          rec(b) = rec(b) + tasa_ion_extraidos(k_e) * prob_migracion(b+3, k_e)
       end do
    end do
    write(*,'(/,A)') ' [N_e VERDADERO] sim = tasa_ion_extraidos vs Fig.3 extraidos:'
    write(*,'(A)')   '   Ne     sim           Fig.3         sim/Fig.3'
    do b = 1, 4
       write(*,'(I5,2ES15.5,F11.3)') b+3, sim(b), pap(b), sim(b)/pap(b)
    end do
    write(*,'(/,A)') ' [N_e RECONSTRUIDO, PE/27] sim = sum_k R_k P(k->j) vs Fig.6:'
    write(*,'(A)')   '   Ne     sim_rec       Fig.6 before   sim/before   sim*eff       Fig.6 after   (sim*eff)/after'
    do b = 1, 4
       write(*,'(I5,2ES15.5,F11.3,2ES15.5,F11.3)') b+3, rec(b), f6_b(b), rec(b)/f6_b(b), &
            rec(b)*eff_ROI(b+3), f6_a(b), rec(b)*eff_ROI(b+3)/f6_a(b)
    end do
    write(*,'(A,F7.4)') '   Retencion global Fig. 6   sum(after)/sum(before)      = ', &
         sum(f6_a) / sum(f6_b)
    write(*,'(A)')      '   (0.2555 = el "75% signal loss in ROI" del texto; cociente de las dos curvas.)'

    ! Volcado para python/roi_cuts_validacion_Xe.py y exporta_datos_tutor.py
    block
      integer :: u_v, u_m
      open(newunit=u_m, file=trim(outdir)//'migracion_Xe.dat', status='replace')
      write(u_m,'(A)') '# j(reconstruido)  k(verdadero)  R_k*P(k->j) [ev/(kg dia)], antes de eff_ROI'
      do b = 1, 4
         do k_e = 1, n_ion
            write(u_m,'(2I5,ES16.7)') b+3, k_e, tasa_ion_extraidos(k_e) * prob_migracion(b+3, k_e)
         end do
      end do
      close(u_m)
      open(newunit=u_v, file=trim(outdir)//'validacion_fig3_fig6_Xe.dat', status='replace')
      write(u_v,'(A)') '# Ne  sim_verdadero  sim_reconstruido  sim_rec_x_effROI  fig3_extraidos  fig6_antes  fig6_despues'
      do b = 1, 4
         write(u_v,'(I3,6ES16.7)') b+3, sim(b), rec(b), rec(b)*eff_ROI(b+3), pap(b), f6_b(b), f6_a(b)
      end do
      close(u_v)
    end block
  end block
  write(*,'(/,A)') ' Espectro en PE:'
  write(*,'(A,F7.1,A,ES13.5)') '   pico en PE = ', (real(ipk,dp)-0.5_dp)*bin_width_pe, ' , valor = ', pe_pk_val
  write(*,'(A,ES13.5)') '   Sum_PE(Total) sobre todo el espectro = ', sum_pe_total
  write(*,'(A,F8.4,A)') '   Sum_PE(Total) / Sum tasa_ion(1:15)   = ', &
       sum_pe_total / max(sum_ext_all, 1.0e-30_dp), '   (~1 si las gaussianas caben en 0-2000 PE)'
  write(*,'(A,/)') '======================================================================================='

  write(*,*) "=== EXITOSO: espectro PE teorico (Total 1SE..15SE), Eventos/(5PE*kg*dia) ==="
end program red100PE_detallado
