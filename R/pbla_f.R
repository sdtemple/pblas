#' f-based PBLA
#'
#' Via 3.3.1., compute pair-based likelihood approximation.
#'
#' @param r numeric vector of removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_f = function(r, beta, gamma, lag = 0){
  if((any(beta <= 0)) | (any(gamma <= 0))){
    return(1e9) # positive rates
  } else{
    # initialize
    n = length(r)
    N = ncol(beta)
    r1 = r[1]
    B = apply(beta[(n+1):N,1:n], 2, sum)
    # calculate log likelihood (line 6)
    ia = rep(-log(n), n)
    ip = - gamma * (r - r1) # check this
    z = ia + ip
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
        if(rj < rk){
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
    for(alpha in 1:n){z[alpha] = z[alpha] + sum(WY[-alpha])}
    z = log(sum(exp(z)))
    # negative log likelihood
    return(-z)
  }
}
