!=======================================================================
! Módulo: mod_stats
! Propósito: Herramientas estadísticas para resolución de detectores.
!=======================================================================
module mod_stats
  use constants, only: dp
  implicit none

contains

  !---------------------------------------------------------------------
  ! Función: poisson_prob
  ! Calcula la probabilidad de observar 'k' eventos cuando el valor
  ! esperado es 'lambda'. Usa logaritmos para estabilidad numérica.
  !---------------------------------------------------------------------
  function poisson_prob(k, lambda) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    real(dp) :: p, log_p, log_fact
    integer :: i
    
    ! Manejo del caso físico donde no se espera ningún electrón
    if (lambda <= 0.0_dp) then
       if (k == 0) then
          p = 1.0_dp
       else
          p = 0.0_dp
       end if
       return
    end if
    
    ! Cálculo del logaritmo del factorial: ln(k!)
    log_fact = 0.0_dp
    do i = 2, k
       log_fact = log_fact + log(real(i, dp))
    end do
    
    ! Fórmula Poisson: ln(P) = k*ln(lambda) - lambda - ln(k!)
    log_p = real(k, dp) * log(lambda) - lambda - log_fact
    
    ! Exponenciamos para recuperar la probabilidad real
    p = exp(log_p)
    
  end function poisson_prob

	! Dentro de mod_stats.f90, añade después de poisson_prob:
  function binomial_prob(k, n, p) result(prob)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: p
    real(dp) :: prob, log_comb, log_p
    integer :: i
    if (p < 0.0_dp .or. p > 1.0_dp .or. k < 0 .or. k > n) then
       prob = 0.0_dp
       return
    end if
    ! Calcular log del coeficiente binomial ln(C(n,k))
    log_comb = 0.0_dp
    do i = 1, k
       log_comb = log_comb + log(real(n - k + i, dp) / real(i, dp))
    end do
    ! log(P) = log_comb + k*log(p) + (n-k)*log(1-p)
    if (k > 0) then
       log_p = log_comb + real(k, dp)*log(p) + real(n-k, dp)*log(1.0_dp - p)
    else
       log_p = real(n, dp)*log(1.0_dp - p)
    end if
    prob = exp(log_p)
  end function binomial_prob

end module mod_stats
