#' Separated Product PBLA (General SEM)
#'
#' Based on product independence, compute pair-based likelihood approximation. Supports exponential infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta numeric rate
#' @param gamma numeric rate
#' @param N integer population size
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_sep_gsem = function(r, beta, gamma, N, A = 1, lag = 0){

  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }

  if((any(beta <= 0)) | (any(gamma <= 0)) |
     (!is.wholenumber(N)) | (N <= 0) |
     (!is.wholenumber(A)) | (A <= 0)){
    # invalid parameters
    return(1e15)
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

    # calculate log likelihood (line 6)
    ia = rep(-log(A), A)
    ip = - delta * (r[1:A] - r1)
    z = ia + ip

    # evaluate psi and chi terms
    XY = rep(0, n)
    b = beta
    denom1 = 2 * delta
    denom2 = b + delta
    for(j in (1:n)){
      X = 0
      Y = 0
      rj= r[j]
      for(k in (1:n)[-j]){
        # lemma 1
        rk = r[k]
        if(rj - lag < rk){
          w = exp(- delta * (rk - rj + lag))
          x = b * delta / denom1 * w
          y = 1 - b * delta / denom1 / denom2 * w
        } else{
          w = exp(- delta * (rj - lag - rk))
          x = b * delta / denom1 * w
          y = delta * (1 + b / denom1 * w) / denom2
        }
        # line twelve
        X = X + x
        Y = Y + log(y)
      }
      XY[j] = log(X * exp(Y))
    }

    # line eight
    for(alpha in 1:A){z[alpha] = z[alpha] + sum(XY[-alpha])}
    z = matrixStats::logSumExp(z)
    a = n * log(gamma / delta)

    # negative log likelihood
    return(-(a+z))
  }
}
