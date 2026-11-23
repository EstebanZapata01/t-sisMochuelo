!=======================================================================
! Programa: chi2_ON_OFF_1D
!   Ajuste chi2 de la amplitud A de CEvNS al residuo ON-OFF de RED-100
!   (arXiv:2411.18641). Como en el articulo, A es el unico parametro:
!     chi2(A) = Sum (dN_i - A R_i)^2/s_i^2 = S3 - 2 A S1 + A^2 S2
!     A_best = S1/S2,   A_90 = A_best + sqrt(2.706/S2)
!
!   Entrada: ionization_spectra_detallado.dat (de red100PE.f90) y
!            red100_residuo_ONOFF_fig8.csv
!   Salida : chi2_ON_OFF_perfil.dat, chi2_ON_OFF_banda.dat,
!            generic_input_Xe.dat, sensib_real_Xe.dat
!=======================================================================
program chi2_ON_OFF_1D
  use constants, only: dp, Z_Ge, N_Ge
  use mod_detector, only: eps_ROI_pe, PE_ROI_min, PE_ROI_max
  implicit none

  integer, parameter :: NMAX = 100

  ! Variables de datos
  integer  :: n_datos
  real(dp) :: pe_dat(NMAX), dN_dat(NMAX), sigma_dat(NMAX)

  ! Prediccion SM (del dat de Fortran)
  integer  :: n_pred
  real(dp), allocatable :: pe_pred(:), R_pred(:)

  ! Prediccion interpolada en bins de datos
  real(dp) :: R_SM(NMAX)

  ! Parametros del analisis
  real(dp), parameter :: dchi2_90 = 2.706_dp
  real(dp), parameter :: exposure_kgd = 192.0_dp   ! exposicion de 2024 en el volumen fiducial
  integer,  parameter :: N_scan   = 10000

  ! Cantidades intermedias
  real(dp) :: S1, S2, S3
  real(dp) :: A
  real(dp) :: chi2_sin
  real(dp) :: chi2_min_sin
  real(dp) :: A_best_sin
  real(dp) :: A_90_sin
  real(dp) :: slope, pe_temp, cols(9)
  real(dp) :: w_dat, bin_pred

  ! Archivos
  integer  :: u_pred, u_perfil, u_banda, u_dat
  real(dp) :: dummy
  integer  :: ios, i, j, j0
  character(len=256) :: line
  character(len=250) :: datadir
  character(len=250) :: f_pred, f_perfil, f_banda

  ! ==================================================================
  ! 0. RUTAS
  ! ==================================================================
  datadir  = '../../datos/'
  f_pred   = trim(datadir)//'ionization_spectra_detallado.dat'
  f_perfil = trim(datadir)//'chi2_ON_OFF_perfil.dat'
  f_banda  = trim(datadir)//'chi2_ON_OFF_banda.dat'

  ! ==================================================================
  ! 1. RESIDUO ON-OFF DE RED-100 (panel E, en PE), leido del PDF del articulo
  ! ==================================================================
  n_datos = 0
  open(newunit=u_dat, file=trim(datadir)//'red100_residuo_ONOFF_fig8.csv', status='old', action='read')
  do
    read(u_dat, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == 'E') then
      n_datos = n_datos + 1
      read(line(3:), *) pe_dat(n_datos), dN_dat(n_datos), sigma_dat(n_datos), dummy
    end if
  end do
  close(u_dat)

  write(*,'(A,I3,A)') '  [1] Cargados ', n_datos, ' bins ON-OFF (extraidos del PDF)'

  ! ==================================================================
  ! 2. LEER PREDICCION SM (espectro en PE de red100PE.f90)
  !    Se usa el total, por la eficiencia del bin reconstruido y por el
  !    cociente de anchos de bin (datos / rejilla).
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
    if (ios /= 0) then
      write(*,*) 'Error al parsear linea: ', trim(line)
      cycle
    end if
    j = j + 1
    pe_pred(j) = cols(1)
    ! cols: 1=PE_center 2=Total 3..9 = 1SE..7SE
    R_pred(j)  = cols(2)          ! Total = suma de las plantillas k=1..15
  end do
  close(u_pred)
  bin_pred = pe_pred(2) - pe_pred(1)                                   ! ancho de la rejilla de red100PE [PE]
  w_dat    = (PE_ROI_max - PE_ROI_min) / real(n_datos, dp)             ! ancho de bin del histograma de datos [PE]
  write(*,'(A,I5,A)') '  [2] Leidos ', n_pred, ' bins de la prediccion SM'

  ! ==================================================================
  ! 3. INTERPOLAR PREDICCION A LOS BINS DE DATOS
  ! ==================================================================
  write(*,'(A)') '  [3] Interpolando prediccion a bins de datos...'
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
      ! Total interpolado (por bin de rejilla) -> por bin de datos, por eff_ROI del bin reconstruido
      R_SM(i) = (R_pred(j0) + slope * (pe_dat(i) - pe_pred(j0))) * (w_dat / bin_pred) &
              * eps_ROI_pe(pe_dat(i))
    else
      write(*,'(A,F7.1,A)') '  AVISO: bin PE=', pe_dat(i), &
        ' fuera del rango de la prediccion -> R_SM=0'
    end if
  end do

  write(*,'(/,A)') '  Verificacion: datos vs prediccion SM interpolada'
  write(*,'(A)')   '  PE_center      dN_data       sigma        R_SM'
  do i = 1, n_datos
    write(*,'(4(F12.4))') pe_dat(i), dN_dat(i), sigma_dat(i), R_SM(i)
  end do

  ! ==================================================================
  ! 3b. ENTRADA PARA chi2_nsi_generic
  ! ==================================================================
  block
    integer :: u_gen, ig
    ! sigma_alpha = -1: sin parametro de molestia
    open(newunit=u_gen, file=trim(datadir)//'generic_input_Xe.dat', status='replace')
    write(u_gen,'(A)') '# Z  N  sigma_alpha  n_bins  ipar'
    write(u_gen,'(2F10.4,F8.4,2I6)') Z_Ge, N_Ge, -1.0_dp, n_datos, 5
    write(u_gen,'(A)') '# bin  dN_dat  sigma_dat  R_SM'
    do ig = 1, n_datos
       write(u_gen,'(I5,3ES16.7)') ig, dN_dat(ig), sigma_dat(ig), R_SM(ig)
    end do
    close(u_gen)
    write(*,'(A)') '  [3b] Volcado generic_input_Xe.dat (validacion cruzada, aditivo)'
  end block

  ! ==================================================================
  ! 4. SUMAS AUXILIARES (independientes de A)
  !    S1 = Sum dNi*Ri/si^2
  !    S2 = Sum Ri^2/si^2
  !    S3 = Sum dNi^2/si^2
  ! ==================================================================
  S1 = 0.0_dp; S2 = 0.0_dp; S3 = 0.0_dp
  do i = 1, n_datos
    S1 = S1 + dN_dat(i) * R_SM(i)  / sigma_dat(i)**2
    S2 = S2 + R_SM(i)**2            / sigma_dat(i)**2
    S3 = S3 + dN_dat(i)**2          / sigma_dat(i)**2
  end do

  ! ==================================================================
  ! 5a. RESULTADO. A_best = max(S1/S2, 0), porque la amplitud no puede
  !     ser negativa; si A_best = 0, A_90 = [S1 + sqrt(S1^2 + 2.706 S2)]/S2
  ! ==================================================================
  block
    real(dp) :: A_best_raw
    A_best_raw = S1 / S2
    A_best_sin = max(A_best_raw, 0.0_dp)
    if (A_best_raw > 0.0_dp) then
      chi2_min_sin = S3 - S1**2 / S2
      A_90_sin     = A_best_sin + sqrt(dchi2_90 / S2)
    else
      chi2_min_sin = S3                                   ! chi2(A=0)
      A_90_sin     = (S1 + sqrt(S1**2 + dchi2_90*S2)) / S2
    end if

    write(*,'(/,A)') '  +--------------------------------------------------+'
    write(*,'(A)')   '  |   RESULTADO (RED-100 SVI, limite observado)      |'
    write(*,'(A)')   '  |   ajuste de 1 solo parametro: la amplitud A      |'
    write(*,'(A)')   '  +--------------------------------------------------+'
    write(*,'(A,F12.5,A)') '  |  A_best (S1/S2, sin restringir) = ', A_best_raw
    write(*,'(A,F12.5,A)') '  |  A_best (fisico, >= 0)          = ', A_best_sin
    write(*,'(A,F12.5,A)') '  |  chi2_min                       = ', chi2_min_sin
    write(*,'(A,F12.5,A)') '  |  A_90% (sup)                    = ', A_90_sin, '   (OBSERVADO)'
    write(*,'(A)')   '  +--------------------------------------------------+'

    if (A_best_raw < 0.0_dp) then
      write(*,'(/,A)') '  -> Mejor ajuste en A<0: senal CEvNS no requerida por los datos'
    else if (A_best_raw <= 1.5_dp) then
      write(*,'(/,A)') '  -> Resultado compatible con la prediccion SM (A~1)'
    end if
  end block

  ! ==================================================================
  ! 5c. SENSIBILIDAD ESPERADA (Asimov, dN_i = R_i): A_90 = 1 + sqrt(2.706/S2),
  !     con las mismas sigma_i del residuo publicado
  ! ==================================================================
  write(*,'(/,A)')     '  +--------------------------------------------------+'
  write(*,'(A)')       '  |   SENSIBILIDAD ESPERADA (Asimov, RED-100 SV)     |'
  write(*,'(A)')       '  +--------------------------------------------------+'
  write(*,'(A,F10.5)') '  |  A_90 esperado = ', 1.0_dp + sqrt(dchi2_90 / S2)
  write(*,'(A,F10.5)') '  |  /sqrt(3) (ajuste a 3 histogramas) = ', &
                          (1.0_dp + sqrt(dchi2_90 / S2)) / sqrt(3.0_dp)
  write(*,'(A)')       '  +--------------------------------------------------+'

  ! ==================================================================
  ! 5b. PERFIL chi2(A)
  ! ==================================================================
  open(newunit=u_perfil, file=f_perfil, status='replace')
  write(u_perfil,'(A)') '# A   chi2(A)'
  do i = 0, N_scan - 1
    A = 0.0_dp + real(i, dp) * (300.0_dp / real(N_scan - 1, dp))
    chi2_sin = S3 - 2.0_dp*A*S1 + A**2*S2
    write(u_perfil,'(2(ES14.6,2X))') A, chi2_sin
  end do
  close(u_perfil)

  ! ==================================================================
  ! VALIDACION: chi2 por bin y compatibilidad con el SM
  ! ==================================================================
  block
    real(dp) :: chi2_1, dchi2_SM, sumRSM, sumAbsdN, contrib
    integer  :: nb
    chi2_1   = S3 - 2.0_dp*S1 + S2          ! chi2(A=1)
    dchi2_SM = chi2_1 - chi2_min_sin
    sumRSM = sum(R_SM(1:n_datos)); sumAbsdN = sum(abs(dN_dat(1:n_datos)))

    write(*,'(/,A)') '  =============== VALIDACION EXTRA ==============='
    write(*,'(A,F9.3,A,I0)') '   chi2_min / ndof        = ', chi2_min_sin/real(n_datos,dp), &
         '   ndof = ', n_datos
    write(*,'(A,3ES12.4)')   '   S1, S2, S3             = ', S1, S2, S3
    write(*,'(A,ES12.4,A)')  '   Sum R_SM (pred. ROI)   = ', sumRSM, ' counts/(kg dia)'
    write(*,'(A,ES12.4)')    '   Sum |dN_dat| (datos)   = ', sumAbsdN
    write(*,'(A,F8.3)')      '   Delta chi2 en SM (A=1) = ', dchi2_SM
    write(*,'(A,F6.2,A)')    '   -> SM compatible a     ~ ', sqrt(max(dchi2_SM,0.0_dp)), ' sigma'
    write(*,'(A)')           '   contribucion por bin a chi2 con A=A_best (buscar outliers de digitalizacion):'
    write(*,'(A)')           '     bin   PE_center   (dN/sigma)^2   R_SM'
    do nb = 1, n_datos
      contrib = ((dN_dat(nb) - A_best_sin*R_SM(nb))/sigma_dat(nb))**2
      write(*,'(I6,F12.3,F14.3,ES13.4,A)') nb, pe_dat(nb), contrib, R_SM(nb), &
           merge(' <-- >4', '       ', contrib > 4.0_dp)
    end do
    write(*,'(A)')           '  ==============================================='
  end block

  ! ==================================================================
  ! 5d. A_90 esperado frente a la exposicion: multiplicar la exposicion por
  !     mult da S2 -> mult*S2, asi que A_90 = 1 + sqrt(2.706/(mult S2))
  ! ==================================================================
  block
    integer, parameter :: n_mult = 8
    real(dp) :: mult_arr(n_mult), A90_esp
    integer  :: im, u_real
    mult_arr = [1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 50.0_dp, 100.0_dp, 335.0_dp, 1000.0_dp]

    open(newunit=u_real, file=trim(datadir)//'sensib_real_Xe.dat', status='replace')
    write(u_real,'(A)') '# mult  exposicion_kgd  A_90_esperado'
    write(*,'(/,A)') '  --- Proyeccion esperada vs exposicion (mult, expo_kgd, A90_esperado) ---'
    do im = 1, n_mult
      A90_esp = 1.0_dp + sqrt(dchi2_90 / (mult_arr(im) * S2))
      write(u_real,'(F8.1,2F16.5)') mult_arr(im), exposure_kgd*mult_arr(im), A90_esp
      write(*,'(F8.1,2F16.5)')      mult_arr(im), exposure_kgd*mult_arr(im), A90_esp
    end do
    close(u_real)
  end block

  ! ==================================================================
  ! 6. GUARDAR LA BANDA A_90*R_SM JUNTO A LOS DATOS
  ! ==================================================================
  open(newunit=u_banda, file=f_banda, status='replace')
  write(u_banda,'(A)') &
    '# PE_center  R_SM  A90*R_SM  -A90*R_SM  delta_ON_OFF  sigma_stat'
  write(u_banda,'(A,F12.4)') '# A_best = ', A_best_sin
  write(u_banda,'(A,F12.4)') '# A_90   = ', A_90_sin
  do i = 1, n_datos
    write(u_banda,'(6(ES14.6,2X))') &
      pe_dat(i), R_SM(i),           &
      A_90_sin * R_SM(i),           &
     -A_90_sin * R_SM(i),           &
      dN_dat(i), sigma_dat(i)
  end do
  close(u_banda)

  write(*,'(/,A)') '=== Archivos generados ==='
  write(*,'(A)') '  chi2_ON_OFF_perfil.dat  -> perfil chi2(A)'
  write(*,'(A)') '  chi2_ON_OFF_banda.dat   -> banda naranja para la figura'

  deallocate(pe_pred, R_pred)

end program chi2_ON_OFF_1D
