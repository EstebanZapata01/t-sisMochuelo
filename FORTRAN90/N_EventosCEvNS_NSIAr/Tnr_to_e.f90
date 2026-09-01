!=======================================================================
! Módulo: mod_tnr_to_e (NEST 218 V/cm - ARGON - lectura robusta y dinámica)
! Rol     : traduce retroceso T [keV] -> electrones via la tabla
!           datos/nest_Ar_218V_dense.txt  (3 col: T, Qy, F), generada con
!           python/nest_Ar.py usando el modelo LArNEST (NR de argon).
! Pipeline: etapa "retroceso -> ionizacion"; misma interfaz que la de Xe
!           (obtener_electrones_creados / obtener_fano / obtener_nest_binomial).
! Tesis   : metodologia.tex Sec. 3.2 y Sec. 9.
! Decision metodologica clave: LArNEST es sub-Poissoniano para NR
!   (F ~ 0.11-0.15 en la ROI, ~0.65 en el umbral); mismo modelo binomial
!   de Fano que Xe. NO se usa F=1.
!=======================================================================
module mod_tnr_to_e
  use constants
  implicit none
  
  ! Arreglos dinámicos que se dimensionarán en tiempo de ejecución
  ! T_dense  : energia de retroceso [keV]
  ! Qy_dense : rendimiento de carga medio [e-/keV]
  ! F_dense  : factor tipo Fano  Var(N_e)/<N_e>  (de GetQuanta de NEST; ~0.43-0.47)
  real(dp), allocatable, dimension(:) :: T_dense, Qy_dense, F_dense
  integer :: n_puntos_actual = 0
  logical :: nest_inicializado = .false.

contains

  !---------------------------------------------------------------------
  subroutine inicializar_nest()
    integer :: u_in, iostat, i, ios3
    character(len=500) :: linea
    real(dp) :: T_temp, Qy_temp, F_temp
    character(len=200) :: archivo_datos

    archivo_datos = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/nest_Ar_218V_dense.txt'

    open(newunit=u_in, file=trim(archivo_datos), status='old', iostat=iostat)
    if (iostat /= 0) then
       print *, "ERROR: No se encontró ", trim(archivo_datos)
       stop
    end if

    ! ---- 1ª pasada: contar líneas de datos (ignorando comentarios) ----
    n_puntos_actual = 0
    do
       read(u_in, '(A)', iostat=iostat) linea
       if (iostat /= 0) exit                ! fin de archivo
       linea = adjustl(linea)
       if (linea(1:1) == '#' .or. len_trim(linea) == 0) cycle  ! ignorar comentarios/vacías
       n_puntos_actual = n_puntos_actual + 1
    end do
    rewind(u_in)

    ! Reservar memoria exacta
    allocate(T_dense(n_puntos_actual), Qy_dense(n_puntos_actual), F_dense(n_puntos_actual))

    ! ---- 2ª pasada: leer los datos (3 col: T Qy F; si solo hay 2, F=1) ----
    i = 0
    do
       read(u_in, '(A)', iostat=iostat) linea
       if (iostat /= 0) exit
       linea = adjustl(linea)
       if (linea(1:1) == '#' .or. len_trim(linea) == 0) cycle
       i = i + 1
       read(linea, *, iostat=ios3) T_temp, Qy_temp, F_temp
       if (ios3 /= 0) then
          read(linea, *) T_temp, Qy_temp
          F_temp = 1.0_dp
       end if
       T_dense(i)  = T_temp
       Qy_dense(i) = Qy_temp
       F_dense(i)  = F_temp
    end do
    close(u_in)

    nest_inicializado = .true.
    print *, "Módulo NEST (218 V/cm): cargados ", n_puntos_actual, " puntos exitosamente."
  end subroutine inicializar_nest

  !---------------------------------------------------------------------
  function obtener_electrones_creados(T_req) result(n_creados)
    real(dp), intent(in) :: T_req
    real(dp) :: n_creados, Qy_calc, pendiente
    integer :: izq, der, medio

    if (.not. nest_inicializado) then
       print *, "ERROR: Debes llamar a inicializar_nest() al principio."
       stop
    end if

    if (T_req <= T_dense(1)) then
       Qy_calc = Qy_dense(1)
    else if (T_req >= T_dense(n_puntos_actual)) then
       Qy_calc = Qy_dense(n_puntos_actual)
    else
       izq = 1; der = n_puntos_actual
       do while (der - izq > 1)
          medio = (izq + der) / 2
          if (T_req >= T_dense(medio)) then
             izq = medio
          else
             der = medio
          end if
       end do
       pendiente = (Qy_dense(der) - Qy_dense(izq)) / (T_dense(der) - T_dense(izq))
       Qy_calc = Qy_dense(izq) + pendiente * (T_req - T_dense(izq))
    end if

    n_creados = T_req * Qy_calc
  end function obtener_electrones_creados

  !---------------------------------------------------------------------
  ! obtener_fano: F(T) = Var(N_e)/<N_e> de NEST, interpolado linealmente.
  !---------------------------------------------------------------------
  function obtener_fano(T_req) result(F_calc)
    real(dp), intent(in) :: T_req
    real(dp) :: F_calc, pendiente
    integer :: izq, der, medio

    if (T_req <= T_dense(1)) then
       F_calc = F_dense(1)
    else if (T_req >= T_dense(n_puntos_actual)) then
       F_calc = F_dense(n_puntos_actual)
    else
       izq = 1; der = n_puntos_actual
       do while (der - izq > 1)
          medio = (izq + der) / 2
          if (T_req >= T_dense(medio)) then
             izq = medio
          else
             der = medio
          end if
       end do
       pendiente = (F_dense(der) - F_dense(izq)) / (T_dense(der) - T_dense(izq))
       F_calc = F_dense(izq) + pendiente * (T_req - T_dense(izq))
    end if
  end function obtener_fano

  !---------------------------------------------------------------------
  ! obtener_nest_binomial: parametros de la binomial de N_e CREADOS que
  ! reproduce (media, varianza) de NEST para el retroceso T_req [keV]:
  !     N_e_creados  ~ Binomial(n_F, p_F)   con  p_F = 1 - F(T),
  !                                              n_F = nint(<N_e>/p_F)
  ! Los EXTRAIDOS (adelgazado binomial con EEE) quedan entonces:
  !     N_e_extraidos ~ Binomial(n_F, p_F*EEE)
  !---------------------------------------------------------------------
  subroutine obtener_nest_binomial(T_req, n_F, p_F)
    real(dp), intent(in)  :: T_req
    integer,  intent(out) :: n_F
    real(dp), intent(out) :: p_F
    real(dp) :: lambda, F

    if (.not. nest_inicializado) then
       print *, "ERROR: Debes llamar a inicializar_nest() al principio."
       stop
    end if

    lambda = obtener_electrones_creados(T_req)   ! <N_e> creados
    F      = obtener_fano(T_req)
    if (F < 0.02_dp) F = 0.02_dp                 ! guardas de seguridad
    if (F > 0.98_dp) F = 0.98_dp
    p_F = 1.0_dp - F

    if (lambda <= 0.0_dp) then
       n_F = 0
    else
       n_F = nint(lambda / p_F)
       if (n_F < 1) n_F = 1
    end if
  end subroutine obtener_nest_binomial

  ! Opcional: liberar memoria al final del programa
  subroutine finalizar_nest()
    if (allocated(T_dense))  deallocate(T_dense)
    if (allocated(Qy_dense)) deallocate(Qy_dense)
    if (allocated(F_dense))  deallocate(F_dense)
    nest_inicializado = .false.
  end subroutine finalizar_nest

end module mod_tnr_to_e
