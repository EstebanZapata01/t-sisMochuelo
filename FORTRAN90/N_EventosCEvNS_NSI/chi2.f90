!=======================================================================
! Programa: chi2_ON_OFF_1D.f90
! Propósito: Leer los datos ON-OFF digitalizados de la Fig. 8 de
!            arXiv:2411.18641 y hacer el análisis chi² 1D sobre la
!            amplitud de señal CEvNS.
!
! Método (RED-100 ajusta solo la amplitud A; sin nuisance de flujo):
!   chi²(A) = Σ (ΔNᵢ - A*Rᵢ)² / σᵢ²  =  S3 - 2*A*S1 + A²*S2
!
!   Donde:
!     A      = amplitud de la señal CEvNS  (A=1 → predicción SM exacta)
!     ΔNᵢ   = residuo ON-OFF digitalizado del artículo [cuentas/kg/día]
!     σᵢ    = incertidumbre estadística en cada bin  [cuentas/kg/día]
!     Rᵢ    = predicción CEvNS SM en bin i  [cuentas/kg/día]
!
!   A_best = S1/S2 ; chi²_min = S3 - S1²/S2 ;
!   A_90   = A_best + sqrt(2.706/S2)   (límite superior, una cola)
!
! Archivos de entrada:
!   ionization_spectra_detallado.dat → predicción de red100PE.f90
!
! Archivos de salida:
!   chi2_ON_OFF_perfil.dat → A  chi2(A)
!   chi2_ON_OFF_banda.dat  → PE  R_SM  banda_sup  banda_inf  datos  sigma
!=======================================================================
program chi2_ON_OFF_1D
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: NMAX = 100   ! máximo de bins de datos

  ! ── Variables de datos ─────────────────────────────────────────────
  integer  :: n_datos
  real(dp) :: pe_dat(NMAX), dN_dat(NMAX), sigma_dat(NMAX)

  ! ── Predicción SM (del dat de Fortran) ────────────────────────────
  integer  :: n_pred
  real(dp), allocatable :: pe_pred(:), R_pred(:)

  ! ── Predicción interpolada en bins de datos ───────────────────────
  real(dp) :: R_SM(NMAX)

  ! ── Parámetros del análisis ────────────────────────────────────────
  real(dp), parameter :: dchi2_90 = 2.706_dp   ! Δchi² → 90% CL (1 cola)
  integer,  parameter :: N_alpha  = 10000      ! puntos del barrido en A

  ! ── Cantidades intermedias ─────────────────────────────────────────
  real(dp) :: S1, S2, S3            ! sumas para la minimización analítica
  real(dp) :: A, dA                 ! amplitud de señal y paso del barrido
  real(dp) :: A_min, A_max
  real(dp) :: chi2_sin
  real(dp) :: chi2_min_sin
  real(dp) :: A_best_sin
  real(dp) :: A_90_sin

  ! ── Archivos ──────────────────────────────────────────────────────
  integer  :: u_pred, u_perfil, u_banda
  integer  :: ios, i, j, j0
  character(len=250) :: datadir
  character(len=250) :: f_pred, f_perfil, f_banda

  ! ── Variables adicionales para lectura robusta ─────────────────────
  character(len=256) :: line
  real(dp) :: pe_temp, R_temp

  ! ══════════════════════════════════════════════════════════════════
  ! 0. RUTAS  (ajusta solo esta línea)
  ! ══════════════════════════════════════════════════════════════════
  datadir = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/'

  f_pred    = trim(datadir) // 'ionization_spectra_detallado.dat'
  f_perfil  = trim(datadir) // 'chi2_ON_OFF_perfil.dat'
  f_banda   = trim(datadir) // 'chi2_ON_OFF_banda.dat'

  ! ══════════════════════════════════════════════════════════════════
  ! 1. DATOS DIGITALIZADOS DEL ARTÍCULO (Ingresados directamente)
  ! ══════════════════════════════════════════════════════════════════
  n_datos = 15
  
  pe_dat(1) = 112.60128_dp; dN_dat(1) = -0.21102_dp; sigma_dat(1) = 0.20630_dp
  pe_dat(2) = 117.86781_dp; dN_dat(2) =  0.09291_dp; sigma_dat(2) = 0.18583_dp
  pe_dat(3) = 123.14301_dp; dN_dat(3) = -0.29449_dp; sigma_dat(3) = 0.18583_dp
  pe_dat(4) = 128.34364_dp; dN_dat(4) =  0.10394_dp; sigma_dat(4) = 0.14961_dp
  pe_dat(5) = 133.57848_dp; dN_dat(5) =  0.14173_dp; sigma_dat(5) = 0.13071_dp
  pe_dat(6) = 138.88895_dp; dN_dat(6) = -0.01732_dp; sigma_dat(6) = 0.12756_dp
  pe_dat(7) = 144.19912_dp; dN_dat(7) = -0.17323_dp; sigma_dat(7) = 0.12913_dp
  pe_dat(8) = 149.46939_dp; dN_dat(8) =  0.09134_dp; sigma_dat(8) = 0.10551_dp
  pe_dat(9) = 154.71783_dp; dN_dat(9) = -0.01417_dp; sigma_dat(9) = 0.09291_dp
  pe_dat(10) = 159.96419_dp; dN_dat(10) = -0.09764_dp; sigma_dat(10) = 0.09291_dp
  pe_dat(11) = 165.24581_dp; dN_dat(11) =  0.04724_dp; sigma_dat(11) = 0.06929_dp
  pe_dat(12) = 170.54956_dp; dN_dat(12) = -0.04094_dp; sigma_dat(12) = 0.07402_dp
  pe_dat(13) = 175.78321_dp; dN_dat(13) =  0.00945_dp; sigma_dat(13) = 0.06614_dp
  pe_dat(14) = 181.02373_dp; dN_dat(14) = -0.01260_dp; sigma_dat(14) = 0.06772_dp
  pe_dat(15) = 186.30357_dp; dN_dat(15) =  0.15118_dp; sigma_dat(15) = 0.05512_dp

  write(*,'(A,I3,A)') '  [1] Cargados ', n_datos, ' bins ON-OFF (Digitalizados)'

  ! ══════════════════════════════════════════════════════════════════
  ! 2. LEER PREDICCIÓN SM (solo columnas 1: PE_center, 2: Total)
  !    Formato esperado del archivo:
  !      PE_center   Total   1SE   2SE   3SE   4SE   5SE   6SE
  !    Se ignoran líneas de comentario (#) y líneas vacías.
  ! ══════════════════════════════════════════════════════════════════
  n_pred = 0
  open(newunit=u_pred, file=f_pred, status='old', action='read')
  ! Primer pase: contar líneas con datos
  do
     read(u_pred, '(A)', iostat=ios) line
     if (ios < 0) exit               ! fin de archivo
     if (ios > 0) cycle              ! error de lectura, saltar
     line = adjustl(line)
     if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
     n_pred = n_pred + 1
  end do
  rewind(u_pred)

  allocate(pe_pred(n_pred), R_pred(n_pred))

  ! Segundo pase: leer los dos primeros números de cada línea
  j = 0
  do
     read(u_pred, '(A)', iostat=ios) line
     if (ios < 0) exit
     if (ios > 0) cycle
     line = adjustl(line)
     if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
     read(line, *, iostat=ios) pe_temp, R_temp
     if (ios /= 0) then
        write(*,*) 'Error al parsear línea: ', trim(line)
        cycle
     end if
     j = j + 1
     pe_pred(j) = pe_temp
     R_pred(j)  = R_temp
  end do
  close(u_pred)
  write(*,'(A,I5,A)') '  [2] Leídos ', n_pred, ' bins de la predicción SM'

  ! ══════════════════════════════════════════════════════════════════
  ! 3. INTERPOLAR PREDICCIÓN A LOS BINS DE DATOS
  !    Interpolación lineal: busca el intervalo [pe_pred(j), pe_pred(j+1)]
  !    que contiene pe_dat(i) y aplica interpolación lineal.
  ! ══════════════════════════════════════════════════════════════════
  write(*,'(A)') '  [3] Interpolando predicción a bins de datos...'
  do i = 1, n_datos
    R_SM(i) = 0.0_dp
    j0 = -1
    do j = 1, n_pred - 1
      if (pe_pred(j) <= pe_dat(i) .and. pe_dat(i) <= pe_pred(j+1)) then
        j0 = j
        exit
      end if
    end do
    if (j0 > 0) then
      ! Interpolación lineal
      R_SM(i) = R_pred(j0) + (R_pred(j0+1) - R_pred(j0)) * (pe_dat(i) - pe_pred(j0)) / (pe_pred(j0+1) - pe_pred(j0))
    else
      write(*,'(A,F7.1,A)') '  ADVERTENCIA: bin PE=', pe_dat(i), &
        ' fuera del rango de la predicción → R_SM=0'
    end if
  end do

  ! Imprimir tabla de verificación
  write(*,'(/,A)') '  Verificación: datos vs predicción SM interpolada'
  write(*,'(A)') '  PE_center      ΔN_data       sigma        R_SM'
  do i = 1, n_datos
    write(*,'(4(F12.4))') pe_dat(i), dN_dat(i), sigma_dat(i), R_SM(i)
  end do

  ! ══════════════════════════════════════════════════════════════════
  ! 4. SUMAS AUXILIARES (independientes de A)
  !    S1 = Σᵢ (ΔNᵢ * Rᵢ) / σᵢ²
  !    S2 = Σᵢ  Rᵢ²        / σᵢ²
  !    S3 = Σᵢ  ΔNᵢ²       / σᵢ²
  ! ══════════════════════════════════════════════════════════════════
  S1 = 0.0_dp; S2 = 0.0_dp; S3 = 0.0_dp
  do i = 1, n_datos
    S1 = S1 + dN_dat(i) * R_SM(i) / sigma_dat(i)**2
    S2 = S2 + R_SM(i)**2           / sigma_dat(i)**2
    S3 = S3 + dN_dat(i)**2         / sigma_dat(i)**2
  end do

  ! ══════════════════════════════════════════════════════════════════
  ! 5a. RESULTADO ANALÍTICO: CHI² SIN NUISANCE DE FLUJO
  !
  !     chi²(A) = Σ (ΔNᵢ - A*Rᵢ)² / σᵢ²
  !             = S3 - 2*A*S1 + A²*S2
  !
  !     Mínimo en A_best = S1/S2
  !     chi²_min = S3 - S1²/S2
  !     Δchi² = S2*(A - A_best)² = 2.706
  !     → A_90 = A_best + sqrt(2.706/S2)   (límite superior una cola)
  ! ══════════════════════════════════════════════════════════════════
  A_best_sin  = S1 / S2
  chi2_min_sin = S3 - S1**2 / S2
  A_90_sin    = A_best_sin + sqrt(dchi2_90 / S2)

  write(*,'(/,A)') '  ┌─────────────────────────────────────────────────┐'
  write(*,'(A)')   '  │   RESULTADO (ajuste de 1 parámetro: A)          │'
  write(*,'(A)')   '  ├─────────────────────────────────────────────────┤'
  write(*,'(A,F10.5,A)') '  │  A_best         = ', A_best_sin,  '                  │'
  write(*,'(A,F10.5,A)') '  │  chi²_min       = ', chi2_min_sin,'                  │'
  write(*,'(A,F10.5,A)') '  │  A_90% (sup)    = ', A_90_sin,    '                  │'
  write(*,'(A)')   '  └─────────────────────────────────────────────────┘'

  if (A_best_sin >= 0.0_dp .and. A_best_sin <= 1.5_dp) then
    write(*,'(/,A)') '  → Resultado compatible con la predicción SM (A~1)'
  else if (A_best_sin < 0.0_dp) then
    write(*,'(/,A)') '  → Mejor ajuste en A<0: señal no requerida por los datos'
  end if

  ! ══════════════════════════════════════════════════════════════════
  ! 5b. PERFIL chi²(A)  (parábola exacta, sin nuisance)
  ! ══════════════════════════════════════════════════════════════════
  A_min = -1.5_dp
  A_max =  5.0_dp
  dA    = (A_max - A_min) / real(N_alpha - 1, dp)

  open(newunit=u_perfil, file=f_perfil, status='replace')
  write(u_perfil,'(A)') '# A   chi2(A)'
  do i = 0, N_alpha - 1
    A = A_min + i * dA
    chi2_sin = S3 - 2.0_dp*A*S1 + A**2*S2
    write(u_perfil,'(2(ES14.6,2X))') A, chi2_sin
  end do
  close(u_perfil)

  ! ══════════════════════════════════════════════════════════════════
  ! 6. GUARDAR BANDA PARA LA FIGURA
  !    Columnas: PE_center  R_SM  banda_sup  banda_inf  datos  sigma
  ! ══════════════════════════════════════════════════════════════════
  open(newunit=u_banda, file=f_banda, status='replace')
  write(u_banda,'(A)') &
    '# PE_center   R_SM   A90*R_SM   -A90*R_SM   delta_ON_OFF   sigma_stat'
  write(u_banda,'(A,F8.5)') '# A_best = ', A_best_sin
  write(u_banda,'(A,F8.5)') '# A_90   = ', A_90_sin
  do i = 1, n_datos
    write(u_banda,'(6(ES14.6,2X))') &
      pe_dat(i), R_SM(i), &
      A_90_sin * R_SM(i), -A_90_sin * R_SM(i), &
      dN_dat(i), sigma_dat(i)
  end do
  close(u_banda)

  write(*,'(/,A)') '=== Archivos generados ==='
  write(*,'(A)') '  chi2_ON_OFF_perfil.dat  → perfil chi²(A), graficar con Python'
  write(*,'(A)') '  chi2_ON_OFF_banda.dat   → datos + banda naranja para la figura'

  deallocate(pe_pred, R_pred)
end program chi2_ON_OFF_1D
