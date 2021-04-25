#' Separated Product PBLA
#'
#' Based on product independence, compute pair-based likelihood approximation. Assume exponential infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_sep = function(r, beta, gamma, lag = 0){
  if((any(beta <= 0)) | (any(gamma <= 0))){
    return(1e9) # positive rates
  } else{
    # initialize
    n = length(r)
    N = ncol(beta)
    r1 = r[1]
    B = apply(beta[(n+1):N,1:n], 2, sum)
    delta = gamma + B
    # calculate log likelihood (line 6)
    ia = rep(-log(n), n)
    ip = - gamma * (r - r1)
    z = ia + ip
    XY = rep(0, n)
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
        if(rj < rk){
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
      XY[j] = log(X * exp(Y))
    }
    # line eight
    for(alpha in 1:n){z[alpha] = z[alpha] + sum(XY[-alpha])}
    z = log(sum(exp(z)))
    a = sum(log(gamma / delta))
    # negative log likelihood
    return(-(a+z))
  }
}
