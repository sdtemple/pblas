#' Product PBLA
#'
#' Based on product independence, compute pair-based likelihood approximation. Assume exponential infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta numeric rate
#' @param gamma numeric rate
#' @param N integer population size
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_prod = function(r, beta, gamma, N, lag = 0){
  if((beta <= 0) | (gamma <= 0)){
    return(1e9)
  } else{
    # initialize
    n = length(r)
    r1 = r[1]
    beta = beta / N
    B = beta * (N - n)
    delta = gamma + B
    # calculate
    ia = rep(-log(n), n)
    ip = - gamma * (r - r1)
    # product expectation (lemma 2)
    pe = sum(log(delta) - log(beta * 1:(n-1) + delta))
    z = ia + ip + pe
    for(j in (1:n)){
      X = 0
      rj = r[j]
      for(k in (1:n)[-j]){
        rk = r[k]
        if(rj < rk){
          x = exp(- delta * (rk - rj + lag))
        } else{
          x = exp(- delta * (rj - lag - rk))
        }
        X = X + x
      }
      z[-j] = z[-j] + log(X) + log(beta) - log(2)
    }
    # line 8
    z = log(sum(exp(z)))
    a = n * (log(gamma) - log(delta))
    # negative log likelihood
    return(-(a+z))
  }
}
