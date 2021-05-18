#' f-based PBLA (General SEM)
#'
#' Via 3.3.1., compute pair-based likelihood approximation. Supports exponential infectious periods.
#'
#' @param r numeric vector of removal times
#' @param beta numeric rate
#' @param gamma numeric rate
#' @param N integer population size
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_f_gsem = function(r, beta, gamma, N, A = 1, lag = 0){

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
    B = beta * (N - n)

    # calculate log likelihood (line 6)
    ia = rep(-log(A), A) # discrete uniform patient zero
    ip = - (gamma + B) * (r[1:A] - r1)
    z = ia + ip

    # evaluate psi, chi, and phi terms
    WY = rep(0, n)
    b = beta
    denom1 = gamma + gamma
    denom2 = gamma + b
    for(j in 1:n){
      W = 0
      Y = 0
      rj = r[j]
      for(k in (1:n)[-j]){
        rk = r[k]
        # f lemmas
        if(rj - lag < rk){
          w = b * gamma / (gamma + gamma + B) * exp(- gamma * (rk - rj + lag))
          x = gamma / denom1 / denom2  * exp(- gamma * (rk - rj + lag))
          y = 1 - b * x
        } else{
          w = b * gamma / (gamma + gamma + B) * exp(- gamma * (rj - lag - rk))
          x = exp(- gamma * (rj - lag - rk)) / denom1
          y = gamma * (1 + b * x) / denom2
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
