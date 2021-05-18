#' f-based PBLA
#'
#' Via 3.3.1., compute pair-based likelihood approximation. Supports exponential infectious periods.
#'
#' @param r numeric vector of removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_f = function(r, beta, gamma, A = 1, lag = 0){

  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }

  if((any(beta <= 0)) | (any(gamma <= 0)) |
     (!is.wholenumber(A)) | (A <= 0)){
    # invalid parameters
    return(1e15)
  } else{

    # initialize
    n = length(r)
    N = ncol(beta)
    r1 = r[1]

    # compute B
    if((n < (N - 1)) & (n > 1)){
      B = apply(beta[1:n,(n+1):N], 1, sum)
    } else{ # handles entire population infected
      if(n == N){B = 0}
      if(n == (N - 1)){B = beta[1:(N-1),N]}
      if(n == 1){B = sum(beta[1,2:N])}
    }

    # calculate log likelihood (line 6)
    ia = rep(-log(A), A) # discrete uniform patient zero
    ip = - (gamma + B)[1:A] * (r[1:A] - r1)
    z = ia + ip

    # evaluate psi, chi, and phi terms
    WY = rep(0, n)
    for(j in (1:n)){
      W = 0
      Y = 0
      rj = r[j]
      gammaj = gamma[j]
      Bj = B[j]
      for(k in (1:n)[-j]){
        b = beta[k,j]
        rk = r[k]
        gammak = gamma[k]
        denom1 = gammaj + gammak
        denom2 = gammak + b
        # f lemmas
        if(rj - lag < rk){
          w = b * gammaj / (gammaj + gammak + Bj) * exp(- gammak * (rk - rj + lag))
          x = gammaj / denom1 / denom2  * exp(- gammak * (rk - rj + lag))
          y = 1 - b * x
        } else{
          w = b * gammaj / (gammaj + gammak + Bj) * exp(- gammaj * (rj - lag - rk))
          x = exp(- gammaj * (rj - lag - rk)) / denom1
          y = gammak * (1 + b * x) / denom2
        }
        W = W + w
        Y = Y + log(y)
      }
      WY[j] = log(W) + Y
    }

    for(alpha in 1:A){z[alpha] = z[alpha] + sum(WY[-alpha])}
    z = matrixStats::logSumExp(z)

    # negative log likelihood
    return(-z)
  }
}
