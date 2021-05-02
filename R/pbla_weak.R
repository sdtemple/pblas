#' Weak limit PBLA
#'
#' Based on weak limit result, compute pair-based likelihood approximation. Assumes exponential infectious periods.
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
pbla_weak = function(r, beta, gamma, N, lag = 0){

  if((beta <= 0) | (gamma <= 0)){
    # invalid parameters
    return(1e9)
  } else{

    # initialize
    n = length(r)
    r1 = r[1]
    beta = beta / N

    # change of variable to delta
    if(n < N){
      B = beta * (N - n)
      delta = gamma + B
    } else{ # handles entire population infected
      if(n == N){delta = gamma}
    }

    # calculate log likelihood (line 8)
    ia = rep(-log(n), n)
    ip = - delta * (r - r1)
    # weak limit (lemma 3)
    wl = - beta / (N * delta) * choose(n, 2) +
      (beta ^ 2) / (12 * ((delta * N) ^ 2) * n * (n - 1) * (4 * n - 5))
    z = ia + ip + wl

    # evaluate chi terms
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
    z = matrixStats::logSumExp(z)
    a = n * (log(gamma) - log(delta))

    # negative log likelihood
    return(-(a+z))
  }
}
