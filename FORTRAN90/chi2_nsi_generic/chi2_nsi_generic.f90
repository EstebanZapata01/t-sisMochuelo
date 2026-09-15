!=======================================================================
! Archivo : chi2_nsi_generic.f90
! Rol     : motor GENERICO chi2(A,alpha) + reduccion NSI a A_amp(eps),
!           independiente del experimento (Z, N y los datos de entrada
!           son parametros leidos de archivo, no constantes de modulo).
! Pipeline: NO forma parte de la cadena de simulacion fisica (no calcula
!           flujo, seccion eficaz ni NEST). Consume la SALIDA de esa
!           cadena (dN/R_exp, sigma, R_pred por bin, ya integrados) que
!           chi2.f90 (Xe) y 2pchi2.f90 (CONUS+) escriben con un bloque
!           write() puramente aditivo (ver Sec. de validacion cruzada,
!           metodologia.tex).
! Tesis   : doc/metodologia.tex, seccion "Validacion cruzada: el mismo
!           codigo con datos de Xe y de CONUS+". Objetivo: demostrar,
!           corriendo UN SOLO binario dos veces, que el motor
!           chi2(A,alpha)+NSI es agnostico al experimento.
! Decision metodologica clave:
!   - El bloque select case(ipar) (15 casos) y el chi2 se copian VERBATIM
!     de chi2red100_nest.f90 y 2pchi2.f90 (ya validados contra el paper de
!     RED-100 y contra la Tabla II de CONUS+). No se rederiva ni se
!     "mejora" nada: es deliberadamente el mismo algoritmo, solo con Z, N
!     y los datos de entrada parametrizados en vez de hardcodeados.
!   - El termino de nuisance de flujo (prior gaussiano con sigma_alpha)
!     se activa solo si sigma_alpha > 0 (caso CONUS+: 0.169). Con
!     sigma_alpha <= 0 (centinela; caso RED-100/Xe) el ajuste es de un
!     solo parametro: chi2(A) = S3 - 2*A*S1 + A^2*S2.
!   - chi2.f90 y 2pchi2.f90 NO se modifican en su logica: solo reciben
!     un bloque write() adicional al final que no altera ningun valor ya
!     calculado ni ningun archivo de salida existente.
! Entradas: archivo de texto (ruta en el 1er argumento de linea de
!           comandos) con el formato:
!             # Z  N  sigma_alpha  n_bins  ipar
!               <Z> <N> <sigma_alpha> <n_bins> <ipar>
!             # bin  dN_o_Rexp   sigma   R_pred
!               1    dN(1)       sigma(1) R_pred(1)
!               ...
!               n_bins ...
! Salidas : <prefijo>_resumen.txt   (Z,N,ipar,A_best,A_90,chi2_min,s_best,banda)
!           <prefijo>_sin2theta.dat (sin2theta_W  A(s)  chi2  dchi2)
!           <prefijo>_perfil.dat    (A  chi2_sinNuisance  chi2_conNuisance;
!                                    mismo formato que chi2_ON_OFF_perfil.dat
!                                    de chi2.f90, para graficarlo con el
!                                    mismo estilo -- fig9_chi2_perfil.pdf --
!                                    pero con datos de otro experimento)
!           <prefijo>_nsi2D.dat     (eps_x eps_y chi2 ; grilla 1000x1000)
!           donde <prefijo> es el 2do argumento de linea de comandos.
!=======================================================================
program chi2_nsi_generic
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  character(len=300) :: f_in, prefix, f_resumen, f_nsi2d
  character(len=256) :: line
  integer :: u_in, u_res, u_out, ios, i, j, j0, n_bins, ipar, n_args

  real(dp) :: Z, N, sigma_alpha
  real(dp), allocatable :: dN(:), sigma(:), R_pred(:)

  real(dp) :: S1, S2, S3, QW_SM
  real(dp) :: q_nsi_ee, q_nsi_emu, q_nsi_etau, q_eff2, A_amp, alpha_best, chi2v
  real(dp) :: eps_x, eps_y, eps_min, eps_max, deps
  real(dp) :: chi2_1d_min, A_1d_min, Aa, cc
  integer, parameter :: n_u = 1000, n_d = 1000

  ! -------------------------------------------------- argumentos / IO
  n_args = command_argument_count()
  if (n_args < 2) then
     write(*,'(A)') 'Uso: chi2_nsi_generic <archivo_entrada> <prefijo_salida>'
     stop 1
  end if
  call get_command_argument(1, f_in)
  call get_command_argument(2, prefix)
  f_resumen = trim(prefix)//'_resumen.txt'
  f_nsi2d   = trim(prefix)//'_nsi2D.dat'

  ! -------------------------------------------------- leer cabecera
  open(newunit=u_in, file=trim(f_in), status='old', action='read')
  call leer_siguiente_dato(u_in, line)
  read(line, *) Z, N, sigma_alpha, n_bins, ipar

  allocate(dN(n_bins), sigma(n_bins), R_pred(n_bins))
  do i = 1, n_bins
     call leer_siguiente_dato(u_in, line)
     read(line, *) j, dN(i), sigma(i), R_pred(i)
  end do
  close(u_in)

  write(*,'(A)')        '===== chi2_nsi_generic: motor generico chi2(A,alpha) + NSI ====='
  write(*,'(A,A)')      ' Archivo de entrada     = ', trim(f_in)
  write(*,'(A,F8.2,A,F8.2)') ' Z, N                   = ', Z, ' , ', N
  write(*,'(A,I4)')      ' n_bins                 = ', n_bins
  write(*,'(A,I4)')      ' ipar                   = ', ipar
  if (sigma_alpha > 0.0_dp) then
     write(*,'(A,F8.4)') ' sigma_alpha (nuisance) = ', sigma_alpha
  else
     write(*,'(A)')      ' sigma_alpha            = (<=0) SIN nuisance de flujo'
  end if

  ! -------------------------------------------------- sumas auxiliares
  S1 = 0.0_dp; S2 = 0.0_dp; S3 = 0.0_dp
  do i = 1, n_bins
     S1 = S1 + dN(i) * R_pred(i)  / sigma(i)**2
     S2 = S2 + R_pred(i)**2       / sigma(i)**2
     S3 = S3 + dN(i)**2           / sigma(i)**2
  end do

  QW_SM = -N/2.0_dp + (1.0_dp - 4.0_dp*0.23857_dp)/2.0_dp * Z

  ! -------------------------------------------------- perfil 1D en A
  ! (identico en forma a chi2red100_nest.f90:155-169; util como chequeo
  !  independiente de A_best/A_90 frente al programa especifico)
  chi2_1d_min = 1.0e30_dp; A_1d_min = 0.0_dp
  Aa = 0.0_dp
  do
     cc = chi2_of_A(Aa)
     if (cc < chi2_1d_min) then; chi2_1d_min = cc; A_1d_min = Aa; end if
     Aa = Aa + 0.02_dp
     if (Aa > 400.0_dp) exit
  end do

  block
     real(dp) :: A90_eff
     A90_eff = -1.0_dp; Aa = A_1d_min
     do
        Aa = Aa + 0.02_dp
        if (Aa > 400.0_dp) exit
        if (chi2_of_A(Aa) - chi2_1d_min >= 2.706_dp) then; A90_eff = Aa; exit; end if
     end do

     open(newunit=u_res, file=trim(f_resumen), status='replace')
     write(u_res,'(A)')          '# resumen chi2_nsi_generic'
     write(u_res,'(A,A)')        '# entrada = ', trim(f_in)
     write(u_res,'(A,F10.3)')    'Z            = ', Z
     write(u_res,'(A,F10.3)')    'N            = ', N
     write(u_res,'(A,I6)')       'ipar         = ', ipar
     write(u_res,'(A,F14.4)')    'QW_SM        = ', QW_SM
     write(u_res,'(A,F14.6)')    'A_best       = ', A_1d_min
     write(u_res,'(A,F14.6)')    'chi2_min     = ', chi2_1d_min
     write(u_res,'(A,F14.6)')    'chi2_min_ndof= ', chi2_1d_min/real(n_bins,dp)
     write(u_res,'(A,F14.4)')    'A_90         = ', A90_eff
     close(u_res)

     write(*,'(A,F12.5)') ' Q_W                    = ', QW_SM
     write(*,'(A,F12.4)') ' A_best (perfil 1D)     = ', A_1d_min
     write(*,'(A,F12.4)') ' chi2_min               = ', chi2_1d_min
     write(*,'(A,F12.4,A,F8.4,A)') ' chi2_min/ndof          = ', &
          chi2_1d_min/real(n_bins,dp), '  (ndof=', real(n_bins,dp), ')'
     write(*,'(A,F12.4)') ' A_90 (Dchi2=2.706)     = ', A90_eff
  end block

  ! -------------------------------------------------- perfil fino chi2(A)
  ! Mismo formato/resolucion que chi2_ON_OFF_perfil.dat (chi2.f90): permite
  ! reusar el mismo estilo de figura (fig9_chi2_perfil.pdf) con datos de
  ! CUALQUIER experimento que se le pase por generic_input_*.dat.
  block
    integer, parameter :: n_scan_p = 10000
    real(dp) :: A_min_p, A_max_p, dA_p, Ap, chi2_sin_p, al_p, chi2_con_p
    integer  :: ip, u_perf
    A_min_p = 0.0_dp; A_max_p = 300.0_dp
    dA_p = (A_max_p - A_min_p) / real(n_scan_p - 1, dp)
    open(newunit=u_perf, file=trim(prefix)//'_perfil.dat', status='replace')
    write(u_perf,'(A)') '# A   chi2_sinNuisance   chi2_conNuisance'
    do ip = 0, n_scan_p - 1
       Ap = A_min_p + ip * dA_p
       chi2_sin_p = S3 - 2.0_dp*Ap*S1 + Ap**2*S2
       if (sigma_alpha > 0.0_dp) then
          al_p = Ap*(S1 - Ap*S2) / (1.0_dp/sigma_alpha**2 + Ap**2*S2)
          chi2_con_p = 0.0_dp
          do j0 = 1, n_bins
             chi2_con_p = chi2_con_p + ((dN(j0) - Ap*(1.0_dp+al_p)*R_pred(j0))/sigma(j0))**2
          end do
          chi2_con_p = chi2_con_p + (al_p/sigma_alpha)**2
       else
          chi2_con_p = chi2_sin_p        ! sin nuisance: ambas columnas iguales
       end if
       write(u_perf,'(3(ES14.6,2X))') Ap, chi2_sin_p, chi2_con_p
    end do
    close(u_perf)
    write(*,'(A,A)') ' Perfil fino escrito en: ', trim(prefix)//'_perfil.dat'
  end block

  ! -------------------------------------------------- barrido en sin^2(theta_W)
  ! Reparametrizacion: sin^2(theta_W) entra en la prediccion CEvNS SOLO via
  !   Q_W(s) = -N/2 + (1-4s)/2 * Z ,
  ! de modo que la amplitud efectiva respecto al SM es
  !   A(s) = [Q_W(s)/Q_W(s0)]^2 ,  s0 = 0.23857 .
  ! El perfil chi2(s) = chi2_of_A(A(s)) reutiliza la misma maquinaria ya
  ! validada (sin rederivar nada). Con datos de CONUS+ debe reproducir el
  ! programa chi2_sin2theta (validacion del metodo).
  block
    integer, parameter :: n_s = 3000
    real(dp) :: s0, smin, smax, s, QWs, As, cs, ds
    real(dp) :: cs_min, s_at_min, s_lo, s_hi
    integer  :: iss, u_s
    s0 = 0.23857_dp
    smin = 0.0_dp; smax = 0.5_dp
    cs_min = 1.0e30_dp; s_at_min = s0
    do iss = 0, n_s
       s   = smin + real(iss,dp)*(smax-smin)/real(n_s,dp)
       QWs = -N/2.0_dp + (1.0_dp - 4.0_dp*s)/2.0_dp * Z
       As  = (QWs / QW_SM)**2
       cs  = chi2_of_A(As)
       if (cs < cs_min) then; cs_min = cs; s_at_min = s; end if
    end do
    s_lo = -1.0_dp; s_hi = -1.0_dp
    open(newunit=u_s, file=trim(prefix)//'_sin2theta.dat', status='replace')
    write(u_s,'(A)') '# sin2theta_W   A(s)   chi2   dchi2'
    do iss = 0, n_s
       s   = smin + real(iss,dp)*(smax-smin)/real(n_s,dp)
       QWs = -N/2.0_dp + (1.0_dp - 4.0_dp*s)/2.0_dp * Z
       As  = (QWs / QW_SM)**2
       cs  = chi2_of_A(As)
       ds  = cs - cs_min
       write(u_s,'(4(ES15.6,2X))') s, As, cs, ds
       if (ds <= 2.706_dp) then
          if (s_lo < 0.0_dp) s_lo = s
          s_hi = s
       end if
    end do
    close(u_s)
    open(newunit=u_res, file=trim(f_resumen), status='old', position='append')
    write(u_res,'(A)')       '# --- barrido sin^2(theta_W) ---'
    write(u_res,'(A,F12.6)') 's_best        = ', s_at_min
    write(u_res,'(A,F12.6,A,F12.6,A)') 's_90CL_band   = [ ', s_lo, ' , ', s_hi, ' ]'
    close(u_res)
    write(*,'(A,F10.5)')         ' s_best (sin2thetaW)    = ', s_at_min
    write(*,'(A,F9.5,A,F9.5,A)') ' banda 90% CL           = [', s_lo, ', ', s_hi, ']'
  end block

  ! -------------------------------------------------- barrido NSI 2D
  ! Bloque copiado VERBATIM de chi2red100_nest.f90:222-267 (mapeo ipar
  ! y chi2 con alpha analitico), solo con Z,N,dN,sigma,R_pred genericos.
  eps_min = -1.0_dp; eps_max = 1.0_dp
  deps = (eps_max - eps_min) / real(n_u - 1, dp)

  open(newunit=u_out, file=trim(f_nsi2d), status='replace')
  write(u_out, '(A)') '# eps_x   eps_y   chi2'
  write(*,'(A)') ' Iniciando barrido NSI generico (1000x1000)...'

  do i = 0, n_u-1
     eps_y = eps_min + i * deps
     do j = 0, n_d-1
        eps_x = eps_min + j * deps
        q_nsi_ee = 0.0_dp; q_nsi_emu = 0.0_dp; q_nsi_etau = 0.0_dp

        select case(ipar)
        case(1);  q_nsi_ee   = (2.0_dp*eps_y + eps_x)*Z + (eps_y + 2.0_dp*eps_x)*N
        case(2);  q_nsi_emu  = (2.0_dp*eps_y + eps_x)*Z + (eps_y + 2.0_dp*eps_x)*N
        case(3);  q_nsi_etau = (2.0_dp*eps_y + eps_x)*Z + (eps_y + 2.0_dp*eps_x)*N
        case(4);  q_nsi_ee = (2.0_dp*eps_x)*Z + (eps_x)*N; q_nsi_emu = (2.0_dp*eps_y)*Z + (eps_y)*N
        case(5);  q_nsi_ee = (eps_x)*Z + (2.0_dp*eps_x)*N; q_nsi_emu = (eps_y)*Z + (2.0_dp*eps_y)*N
        case(6);  q_nsi_ee = (2.0_dp*eps_x)*Z + (eps_x)*N; q_nsi_emu = (eps_y)*Z + (2.0_dp*eps_y)*N
        case(7);  q_nsi_ee = (eps_x)*Z + (2.0_dp*eps_x)*N; q_nsi_emu = (2.0_dp*eps_y)*Z + (eps_y)*N
        case(8);  q_nsi_ee = (2.0_dp*eps_x)*Z + (eps_x)*N; q_nsi_etau = (2.0_dp*eps_y)*Z + (eps_y)*N
        case(9);  q_nsi_ee = (eps_x)*Z + (2.0_dp*eps_x)*N; q_nsi_etau = (eps_y)*Z + (2.0_dp*eps_y)*N
        case(10); q_nsi_ee = (2.0_dp*eps_x)*Z + (eps_x)*N; q_nsi_etau = (eps_y)*Z + (2.0_dp*eps_y)*N
        case(11); q_nsi_ee = (eps_x)*Z + (2.0_dp*eps_x)*N; q_nsi_etau = (2.0_dp*eps_y)*Z + (eps_y)*N
        case(12); q_nsi_emu = (2.0_dp*eps_x)*Z + (eps_x)*N; q_nsi_etau = (2.0_dp*eps_y)*Z + (eps_y)*N
        case(13); q_nsi_emu = (eps_x)*Z + (2.0_dp*eps_x)*N; q_nsi_etau = (eps_y)*Z + (2.0_dp*eps_y)*N
        case(14); q_nsi_emu = (2.0_dp*eps_x)*Z + (eps_x)*N; q_nsi_etau = (eps_y)*Z + (2.0_dp*eps_y)*N
        case(15); q_nsi_emu = (eps_x)*Z + (2.0_dp*eps_x)*N; q_nsi_etau = (2.0_dp*eps_y)*Z + (eps_y)*N
        end select

        q_eff2 = (QW_SM + q_nsi_ee)**2 + q_nsi_emu**2 + q_nsi_etau**2
        A_amp  = q_eff2 / QW_SM**2

        if (sigma_alpha > 0.0_dp) then
           alpha_best = A_amp*(S1 - A_amp*S2) / (1.0_dp/sigma_alpha**2 + A_amp**2*S2)
        else
           alpha_best = 0.0_dp
        end if

        chi2v = 0.0_dp
        do j0 = 1, n_bins
           chi2v = chi2v + ((dN(j0) - A_amp*(1.0_dp + alpha_best)*R_pred(j0)) / sigma(j0))**2
        end do
        if (sigma_alpha > 0.0_dp) chi2v = chi2v + (alpha_best/sigma_alpha)**2

        write(u_out, '(3ES15.6)') eps_x, eps_y, chi2v
     end do
     write(u_out, *)
  end do
  close(u_out)

  write(*,'(A,A)') ' Resumen escrito en   : ', trim(f_resumen)
  write(*,'(A,A)') ' Grilla NSI escrita en: ', trim(f_nsi2d)
  write(*,'(A,/)') '=================================================================='

  deallocate(dN, sigma, R_pred)

contains

  function chi2_of_A(A) result(c)
    real(dp), intent(in) :: A
    real(dp) :: c, al
    integer  :: m
    if (sigma_alpha > 0.0_dp) then
       al = A*(S1 - A*S2) / (1.0_dp/sigma_alpha**2 + A**2*S2)
    else
       al = 0.0_dp
    end if
    c = 0.0_dp
    do m = 1, n_bins
       c = c + ((dN(m) - A*(1.0_dp + al)*R_pred(m)) / sigma(m))**2
    end do
    if (sigma_alpha > 0.0_dp) c = c + (al/sigma_alpha)**2
  end function chi2_of_A

  ! salta lineas vacias/comentario ('#') y devuelve la siguiente linea util
  subroutine leer_siguiente_dato(unit, out_line)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: out_line
    integer :: iost
    do
       read(unit, '(A)', iostat=iost) out_line
       if (iost /= 0) stop 'chi2_nsi_generic: fin de archivo inesperado leyendo entrada'
       out_line = adjustl(out_line)
       if (len_trim(out_line) == 0) cycle
       if (out_line(1:1) == '#') cycle
       exit
    end do
  end subroutine leer_siguiente_dato

end program chi2_nsi_generic
