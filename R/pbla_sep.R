#' Separated Product PBLA
#'
#' Based on product independence, compute pair-based likelihood approximation. Supports exponential infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_sep = function(r, beta, gamma, A = 1, lag = 0){

  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }

  if((any(beta < 0)) | (any(gamma < 0)) |
     (!is.wholenumber(A)) | (A <= 0)){
    # invalid parameters
    return(1e15)
  } else{

    # initialize
    n = length(r)
    N = ncol(beta)
    r1 = r[1]

    # change of variable to delta
    if((n < (N - 1)) & (n > 1)){
      B = apply(beta[1:n,(n+1):N], 1, sum)
      delta = gamma + B
    } else{ # handles special cases
      if(n == N){delta = gamma}
      if(n == (N - 1)){delta = gamma + beta[1:(N-1),N]}
      if(n == 1){delta = gamma + sum(beta[1,2:N])}
    }

    # calculate log likelihood (line 6)
    ia = rep(-log(A), A)
    ip = - delta[1:A] * (r[1:A] - r1)
    z = ia + ip

    # evaluate psi and chi terms
    psichi = rep(0, n)
    for(j in (1:n)){
      X = 0
      Y = 0
      rj= r[j]
      deltaj = delta[j]
      for(k in (1:n)[-j]){
        # lemma 1
        b = beta[k,j]
        rk = r[k]
        deltak = delta[k]
        denom1 = (deltaj + deltak)
        denom2 = (b + deltak)
        if(rj - lag < rk){
          w = exp(- deltak * (rk - rj + lag))
          x = b * deltaj / denom1 * w
          y = 1 - b * deltaj / denom1 / denom2 * w
        } else{
          w = exp(- deltaj * (rj - lag - rk))
          x = b * deltaj / denom1 * w
          y = deltak * (1 + b / denom1 * w) / denom2
        }
        # line twelve
        X = X + x
        Y = Y + log(y)
      }
      psichi[j] = log(X * exp(Y))
    }

    # line eight
    for(alpha in 1:A){z[alpha] = z[alpha] + sum(psichi[-alpha])}
    z = matrixStats::logSumExp(z)
    a = sum(log(gamma / delta))

    # negative log likelihood
    return(-(a+z))
  }
}
