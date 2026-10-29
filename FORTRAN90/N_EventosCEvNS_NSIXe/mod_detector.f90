!=======================================================================
! Archivo : mod_detector.f90
! Rol     : respuesta del detector en fotoelectrones (PE): ganancia
!           SEG = 27 PE/e-, resolucion de un electron unico sig1, y
!           conversion N_e -> PE.
! Pipeline: etapa "extraccion / PE" -> la usan red100PE y mainred100_nest
!           para el espectro en energia corregida [PE] que ajusta chi2.f90.
! Tesis   : metodologia.tex Sec. 3.3 y Sec. 4 (ajuste 1D ON-OFF).
! Decision metodologica clave: SEG/sig1 son de LXe; en la carpeta de Ar
!   este modulo NO es fisico (solo diagnostico), ver Sec. 9.
!=======================================================================
module mod_detector
  use constants, only: dp, eff_ROI
  implicit none

  ! ==================== PARÁMETROS FÍSICOS ====================
  real(dp), parameter :: SEG_val = 27.0_dp  ! PE/e-, paper calibración
  real(dp), parameter :: sig1    = 7.6_dp   ! PE, sigma del 1SE (Fig.12)
  real(dp), parameter :: pi_val  = 3.141592653589793_dp

  ! ROI del analisis en energia corregida [PE]: 4 a 7 e- = 110 a 189 PE (paper, Sec. IV)
  real(dp), parameter :: PE_ROI_min = 110.0_dp, PE_ROI_max = 189.0_dp

  integer, parameter :: max_curvas = 7   ! plantillas 1SE..7SE (la ROI del paper llega a N_e=7)

  ! ==================== GRILLA PRE-CALCULADA ====================
  integer,  parameter :: N_GRID  = 5000
  real(dp), parameter :: PE_MIN  = 0.0_dp
  real(dp), parameter :: PE_MAX  = 300.0_dp

  real(dp) :: x_grid(N_GRID)
  real(dp) :: y_gauss(N_GRID, max_curvas)
  logical  :: SE_loaded = .false.

contains

  !-------------------------------------------------------------------
  ! cargar_SE_data: mantiene el mismo nombre y firma que el original.
  ! El argumento 'filename' se ignora — las gaussianas se calculan
  ! analíticamente. El resto del código no necesita cambiar nada.
  !-------------------------------------------------------------------
  subroutine cargar_SE_data(filename)
    character(len=*), intent(in) :: filename
    integer  :: i, k
    real(dp) :: dx, mu, sig, norm_check

    write(*,*) "----------------------------------------------------------"
    write(*,*) "   MOD_DETECTOR: Gaussianas analiticas (arXiv:2411.18641) "
    write(*,*) "   SEG  = ", SEG_val, " PE/e-"
    write(*,*) "   sig1 = ", sig1,    " PE"
    write(*,*) "   (archivo ignorado: ", trim(filename), ")"
    write(*,*) "----------------------------------------------------------"

    dx = (PE_MAX - PE_MIN) / real(N_GRID - 1, dp)
    do i = 1, N_GRID
       x_grid(i) = PE_MIN + (i - 1) * dx
    end do

    do k = 1, max_curvas
       mu  = real(k, dp) * SEG_val
       sig = sqrt(real(k, dp)) * sig1

       do i = 1, N_GRID
          y_gauss(i, k) = exp(-0.5_dp * ((x_grid(i) - mu) / sig)**2) &
                          / (sqrt(2.0_dp * pi_val) * sig)
       end do

       norm_check = sum(y_gauss(:, k)) * dx
       write(*,'(A,I2,A,F6.1,A,F5.2,A,F8.5)') &
            "   k=", k, "SE: mu=", mu, " sig=", sig, &
            "  integral=", norm_check
    end do

    SE_loaded = .true.
    write(*,*) "----------------------------------------------------------"
  end subroutine cargar_SE_data

  !-------------------------------------------------------------------
  ! respuesta_empirica: mismo nombre, misma firma, mismo resultado [PE^-1].
  ! Ahora interpola sobre la gaussiana pre-calculada en lugar del PCHIP.
  !-------------------------------------------------------------------
  function respuesta_empirica(S, k) result(val)
    real(dp), intent(in) :: S
    integer,  intent(in) :: k
    real(dp) :: val
    integer  :: lo, hi, mid

    val = 0.0_dp
    if (.not. SE_loaded .or. k < 1 .or. k > max_curvas) return

    if (S <= x_grid(1))      then; val = y_gauss(1,      k); return; end if
    if (S >= x_grid(N_GRID)) then; val = y_gauss(N_GRID, k); return; end if

    ! Búsqueda binaria (idéntica al original)
    lo = 1; hi = N_GRID
    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       if (S >= x_grid(mid)) then; lo = mid; else; hi = mid; end if
    end do

    ! Interpolación lineal (idéntica al original)
    val = y_gauss(lo, k) + (y_gauss(hi, k) - y_gauss(lo, k)) &
          * (S - x_grid(lo)) / (x_grid(hi) - x_grid(lo))

  end function respuesta_empirica

  !-------------------------------------------------------------------
  ! N_e verdadero k -> PE ~ Normal(SEG*k, sqrt(k)*sig1)  [PE^-1], analitica
  !-------------------------------------------------------------------
  function respuesta_gauss(S, k) result(val)
    real(dp), intent(in) :: S
    integer,  intent(in) :: k
    real(dp) :: val, sig
    sig = sqrt(real(k, dp)) * sig1
    val = exp(-0.5_dp * ((S - real(k, dp) * SEG_val) / sig)**2) / (sqrt(2.0_dp * pi_val) * sig)
  end function respuesta_gauss

  !-------------------------------------------------------------------
  ! N_e RECONSTRUIDO = PE corregido / SEG, al entero mas cercano, dentro de la
  ! ROI (110-189 PE); 0 fuera. Es el N_e del eje de las senales del paper.
  !-------------------------------------------------------------------
  function bin_ne_rec(S) result(n)
    real(dp), intent(in) :: S
    integer :: n
    n = 0
    if (S < PE_ROI_min .or. S > PE_ROI_max) return
    n = max(4, min(7, nint(S / SEG_val)))
  end function bin_ne_rec

  !-------------------------------------------------------------------
  ! eff_ROI del bin reconstruido al que cae S (0 fuera de la ROI)
  !-------------------------------------------------------------------
  function eps_ROI_pe(S) result(eps)
    real(dp), intent(in) :: S
    real(dp) :: eps
    integer  :: n
    n = bin_ne_rec(S)
    eps = 0.0_dp
    if (n >= 4) eps = eff_ROI(n)
  end function eps_ROI_pe

  !-------------------------------------------------------------------
  ! Matriz de migracion: P(N_e verdadero k -> N_e reconstruido j), j = 4..7.
  ! Ventana del bin j: [(j-1/2)*SEG, (j+1/2)*SEG] recortada a la ROI.
  !-------------------------------------------------------------------
  function prob_migracion(j, k) result(p)
    integer, intent(in) :: j, k
    real(dp) :: p, a, b, mu, sg, r2
    r2 = sqrt(2.0_dp)
    a  = max((real(j, dp) - 0.5_dp) * SEG_val, PE_ROI_min)
    b  = min((real(j, dp) + 0.5_dp) * SEG_val, PE_ROI_max)
    mu = real(k, dp) * SEG_val
    sg = sqrt(real(k, dp)) * sig1
    p  = 0.5_dp * (erf((b - mu) / (sg * r2)) - erf((a - mu) / (sg * r2)))
  end function prob_migracion

end module mod_detector
