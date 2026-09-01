!=======================================================================
! Programa: chi2_ON_OFF_1D.f90
! Proposito: Leer los datos ON-OFF digitalizados de arXiv:2403.12645
!            y hacer el analisis chi2 1D sobre la amplitud de senal CEvNS.
!
! Metodo:
!   chi2(A, alpha) = Sum (dNi - A*(1+alpha)*Ri)^2/si^2 + alpha^2/sigF^2
!
!   Donde:
!     A      = amplitud de la senal CEvNS  (A=1 -> prediccion SM exacta)
!     alpha  = nuisance de normalizacion de flujo  (pull term)
!     dNi    = residuo ON-OFF digitalizado del articulo [cuentas/kg/dia]
!     si     = incertidumbre estadistica en cada bin  [cuentas/kg/dia]
!     Ri     = prediccion CEvNS SM en bin i  [cuentas/kg/dia]
!     sigF   = incertidumbre sistematica del flujo = 16.9%
!
!   Para cada valor de A, se minimiza analiticamante sobre alpha.
!   El limite al 90% C.L. es donde Dchi2 = chi2(A) - chi2_min = 2.706
!
! Archivos de entrada:
!   ionization_spectra_detallado.dat  ->  prediccion de red100PE.f90
!
! Archivos de salida:
!   chi2_ON_OFF_perfil.dat  ->  A  chi2_sinNuisance  chi2_conNuisance
!   chi2_ON_OFF_banda.dat   ->  PE R_SM banda_sup banda_inf datos sigma
!
! Pipeline: unica rama con DATOS REALES de RED-100. NO se clona a Ar.
! Tesis   : metodologia.tex Sec. 4 (reproduccion del analisis de Xe, P1).
!           Resultado A_90 ~ 111 xSM; la diferencia con la Tabla I del
!           paper = 1/3 histogramas (sqrt3) + modelo de espectro +
!           aceptancia del detector (Fig.3 vs Fig.6).
! Decision metodologica clave: alpha minimizado analiticamente,
!   alpha_hat(A) = A(S1 - A S2)/(1/sigF^2 + A^2 S2); sin nuisance,
!   A_90 = S1/S2 + sqrt(2.706/S2)  (Ecs. 20-22 de metodologia.tex).
!=======================================================================
program chi2_ON_OFF_1D
  use constants, only: dp, eff_ROI
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
  real(dp), parameter :: sigma_F  = 0.169_dp
  real(dp), parameter :: dchi2_90 = 2.706_dp
  integer,  parameter :: N_scan   = 10000

  ! Cantidades intermedias
  real(dp) :: S1, S2, S3
  real(dp) :: A, dA, A_min, A_max
  real(dp) :: alpha_nu
  real(dp) :: chi2_sin, chi2_con
  real(dp) :: chi2_min_sin, chi2_min_con
  real(dp) :: A_best_sin, A_best_con
  real(dp) :: A_90_sin,   A_90_con
  real(dp) :: slope, pe_temp, cols(9)

  ! Archivos
  integer  :: u_pred, u_perfil, u_banda
  integer  :: ios, i, j, j0
  character(len=256) :: line
  character(len=250) :: datadir
  character(len=250) :: f_pred, f_perfil, f_banda

  ! ==================================================================
  ! 0. RUTAS
  ! ==================================================================
  datadir  = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'
  f_pred   = trim(datadir)//'ionization_spectra_detallado.dat'
  f_perfil = trim(datadir)//'chi2_ON_OFF_perfil.dat'
  f_banda  = trim(datadir)//'chi2_ON_OFF_banda.dat'

  ! ==================================================================
  ! 1. DATOS DIGITALIZADOS DEL ARTICULO (hardcoded)
  ! ==================================================================
  n_datos = 15

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

  write(*,'(A,I3,A)') '  [1] Cargados ', n_datos, ' bins ON-OFF (digitalizados)'

  ! ==================================================================
  ! 2. LEER PREDICCION SM (espectro teorico de red100PE.f90)
  !    Formato (9 col): PE_center  Total  1SE 2SE 3SE 4SE 5SE 6SE 7SE
  !    La prediccion para el chi2 es la suma de las plantillas de la ROI
  !    (N_e = 4..7) pesadas por la eficiencia de cortes eff_ROI(N_e)
  !    (digitalizada de la Fig. 6 de arXiv:2411.18641). Las plantillas
  !    N_e < 4 y N_e > 7 quedan fuera de la ROI.
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
    R_pred(j)  = cols(6)*eff_ROI(4) + cols(7)*eff_ROI(5) &
               + cols(8)*eff_ROI(6) + cols(9)*eff_ROI(7)
  end do
  close(u_pred)
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
      R_SM(i) = R_pred(j0) + slope * (pe_dat(i) - pe_pred(j0))
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
  ! 5a. CHI2 SIN NUISANCE (analitico)
  !     chi2(A) = S3 - 2*A*S1 + A^2*S2
  !     A_best  = S1/S2
  !     chi2_min = S3 - S1^2/S2
  !     A_90 = A_best + sqrt(2.706/S2)
  ! ==================================================================
  A_best_sin   = S1 / S2
  chi2_min_sin = S3 - S1**2 / S2
  A_90_sin     = A_best_sin + sqrt(dchi2_90 / S2)

  write(*,'(/,A)') '  +--------------------------------------------------+'
  write(*,'(A)')   '  |       RESULTADOS SIN NUISANCE DE FLUJO           |'
  write(*,'(A)')   '  +--------------------------------------------------+'
  write(*,'(A,F10.5)') '  |  A_best      = ', A_best_sin
  write(*,'(A,F10.5)') '  |  chi2_min    = ', chi2_min_sin
  write(*,'(A,F10.5,A)') '  |  A_90% (sup) = ', A_90_sin, '   (OBSERVADO, RED-100 SVI)'
  write(*,'(A)')   '  +--------------------------------------------------+'

  ! ==================================================================
  ! 5c. SENSIBILIDAD ESPERADA  (RED-100 SV)
  !   Asimov: el "dato" es la senal SM  =>  dN_i = R_i, A_best = 1.
  !   chi2(A) = (1-A)^2 * sum_i R_i^2/sigma_i^2 = (1-A)^2 * S2
  !   A_90_esperado = 1 + sqrt(2.706/S2).
  !   Usa las MISMAS sigma_i del residuo ON-OFF publicado (que ya llevan la
  !   estadistica del fondo); no se simula ningun fondo para Xe.
  !   Reproduce (salvo el factor sqrt(3) por ajustar 1 histograma en vez de
  !   los 3 de RED-100) la columna "sensibilidad" de la Tabla I de
  !   arXiv:2411.18641 (SM2018 58, KI 90, DB 56, INR 64 xSM). Valida el
  !   metodo SV que chi2_bkg_nest.f90 (Ar) aplica con fondo simulado.
  ! ==================================================================
  write(*,'(/,A)')     '  +--------------------------------------------------+'
  write(*,'(A)')       '  |   SENSIBILIDAD ESPERADA (Asimov, RED-100 SV)     |'
  write(*,'(A)')       '  +--------------------------------------------------+'
  write(*,'(A,F10.5)') '  |  A_90 esperado = ', 1.0_dp + sqrt(dchi2_90 / S2)
  write(*,'(A,F10.5)') '  |  /sqrt(3) (ajuste a 3 histogramas) = ', &
                          (1.0_dp + sqrt(dchi2_90 / S2)) / sqrt(3.0_dp)
  write(*,'(A)')       '  +--------------------------------------------------+'

  ! ==================================================================
  ! 5b. CHI2 CON NUISANCE DE FLUJO (barrido numerico + min analitico)
  !
  !     Minimizando d(chi2)/d(alpha) = 0:
  !       alpha_min(A) = A*(S1 - A*S2) / (1/sigF^2 + A^2*S2)
  ! ==================================================================
  A_min = 0.0_dp
  A_max =  300.0_dp
  dA    = (A_max - A_min) / real(N_scan - 1, dp)

  chi2_min_con = 1.0d30
  A_best_con   = 0.0_dp
  A_90_con     = A_max

  open(newunit=u_perfil, file=f_perfil, status='replace')
  write(u_perfil,'(A)') '# A   chi2_sinNuisance   chi2_conNuisance'

  do i = 0, N_scan - 1
    A = A_min + i * dA

    chi2_sin = S3 - 2.0_dp*A*S1 + A**2*S2

    alpha_nu = A*(S1 - A*S2) / (1.0_dp/sigma_F**2 + A**2*S2)

    chi2_con = 0.0_dp
    do j = 1, n_datos
      chi2_con = chi2_con + &
        ((dN_dat(j) - A*(1.0_dp + alpha_nu)*R_SM(j)) / sigma_dat(j))**2
    end do
    chi2_con = chi2_con + (alpha_nu/sigma_F)**2

    write(u_perfil,'(3(ES14.6,2X))') A, chi2_sin, chi2_con

    if (chi2_con < chi2_min_con) then
      chi2_min_con = chi2_con
      A_best_con   = A
    end if
  end do
  close(u_perfil)

  ! Buscar A_90 con nuisance (cruce Dchi2 = 2.706 hacia arriba)
  do i = 0, N_scan - 1
    A = A_min + i * dA
    if (A < A_best_con) cycle
    alpha_nu = A*(S1 - A*S2) / (1.0_dp/sigma_F**2 + A**2*S2)
    chi2_con = 0.0_dp
    do j = 1, n_datos
      chi2_con = chi2_con + &
        ((dN_dat(j) - A*(1.0_dp + alpha_nu)*R_SM(j))/sigma_dat(j))**2
    end do
    chi2_con = chi2_con + (alpha_nu/sigma_F)**2
    if (chi2_con - chi2_min_con >= dchi2_90) then
      A_90_con = A; exit
    end if
  end do

  write(*,'(/,A)') '  +--------------------------------------------------+'
  write(*,'(A)')   '  |     RESULTADOS CON NUISANCE DE FLUJO 16.9%       |'
  write(*,'(A)')   '  +--------------------------------------------------+'
  write(*,'(A,F10.5)') '  |  A_best      = ', A_best_con
  write(*,'(A,F10.5)') '  |  chi2_min    = ', chi2_min_con
  write(*,'(A,F10.5)') '  |  A_90% (sup) = ', A_90_con
  write(*,'(A)')   '  +--------------------------------------------------+'

  if (A_best_con >= 0.0_dp .and. A_best_con <= 1.5_dp) then
    write(*,'(/,A)') '  -> Resultado compatible con la prediccion SM (A~1)'
  else if (A_best_con < 0.0_dp) then
    write(*,'(/,A)') '  -> Mejor ajuste en A<0: senal no requerida por los datos'
  end if

  ! ==================================================================
  !            VALIDACION EXTRA (vs arXiv:2411.18641)
  ! ==================================================================
  block
    real(dp) :: alpha1, chi2_1, dchi2_SM, sumRSM, sumAbsdN, contrib
    integer  :: nb
    alpha1 = 1.0_dp*(S1 - 1.0_dp*S2) / (1.0_dp/sigma_F**2 + 1.0_dp*S2)
    chi2_1 = 0.0_dp
    do nb = 1, n_datos
      chi2_1 = chi2_1 + ((dN_dat(nb) - 1.0_dp*(1.0_dp+alpha1)*R_SM(nb))/sigma_dat(nb))**2
    end do
    chi2_1 = chi2_1 + (alpha1/sigma_F)**2
    dchi2_SM = chi2_1 - chi2_min_con
    sumRSM = sum(R_SM(1:n_datos)); sumAbsdN = sum(abs(dN_dat(1:n_datos)))

    write(*,'(/,A)') '  =============== VALIDACION EXTRA ==============='
    write(*,'(A,F9.3,A,I0)') '   chi2_min / ndof        = ', chi2_min_con/real(n_datos,dp), &
         '   ndof = ', n_datos
    write(*,'(A,3ES12.4)')   '   S1, S2, S3             = ', S1, S2, S3
    write(*,'(A,ES12.4,A)')  '   Sum R_SM (pred. ROI)   = ', sumRSM, ' counts/(kg dia)'
    write(*,'(A,ES12.4)')    '   Sum |dN_dat| (datos)   = ', sumAbsdN
    write(*,'(A,F8.3)')      '   alpha(A=1)             = ', alpha1
    write(*,'(A,F8.3)')      '   Delta chi2 en SM (A=1) = ', dchi2_SM
    write(*,'(A,F6.2,A)')    '   -> SM compatible a     ~ ', sqrt(max(dchi2_SM,0.0_dp)), ' sigma'
    write(*,'(A)')           '   contribucion por bin a chi2 con A=A_best (buscar outliers de digitalizacion):'
    write(*,'(A)')           '     bin   PE_center   (dN/sigma)^2   R_SM'
    do nb = 1, n_datos
      contrib = ((dN_dat(nb) - A_best_con*(1.0_dp+0.0_dp)*R_SM(nb))/sigma_dat(nb))**2
      write(*,'(I6,F12.3,F14.3,ES13.4,A)') nb, pe_dat(nb), contrib, R_SM(nb), &
           merge(' <-- >4', '       ', contrib > 4.0_dp)
    end do
    write(*,'(A)')           '  ==============================================='
  end block

  ! ==================================================================
  ! 6. GUARDAR BANDA PARA LA FIGURA
  !    Columnas: PE_center R_SM A90*R_SM -A90*R_SM delta_ON_OFF sigma
  ! ==================================================================
  open(newunit=u_banda, file=f_banda, status='replace')
  write(u_banda,'(A)') &
    '# PE_center  R_SM  A90*R_SM  -A90*R_SM  delta_ON_OFF  sigma_stat'
  write(u_banda,'(A,F12.4)') '# A_best (sin nuisance) = ', A_best_sin
  write(u_banda,'(A,F12.4)') '# A_90   (sin nuisance) = ', A_90_sin
  write(u_banda,'(A,F12.4)') '# A_best (con nuisance) = ', A_best_con
  write(u_banda,'(A,F12.4)') '# A_90   (con nuisance) = ', A_90_con
  do i = 1, n_datos
    write(u_banda,'(6(ES14.6,2X))') &
      pe_dat(i), R_SM(i),           &
      A_90_con * R_SM(i),           &
     -A_90_con * R_SM(i),           &
      dN_dat(i), sigma_dat(i)
  end do
  close(u_banda)

  write(*,'(/,A)') '=== Archivos generados ==='
  write(*,'(A)') '  chi2_ON_OFF_perfil.dat  -> perfil chi2(A)'
  write(*,'(A)') '  chi2_ON_OFF_banda.dat   -> banda naranja para la figura'

  deallocate(pe_pred, R_pred)
end program chi2_ON_OFF_1D
