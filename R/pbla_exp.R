#' Standard PBLA (Exponential)
#'
#' Compute pair-based likelihood approximation. Assume exponential infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_exp = function(r, beta, gamma, lag = 0){
  if((any(beta <= 0)) | (any(gamma <= 0))){
    return(1e9) # positive rates
  } else{
    # initialize
    n = length(r)
    N = ncol(beta)
    r1 = r[1]
    B = apply(beta[(n+1):N,1:n], 2, sum)
    delta = gamma + B
    # calculate expectations (line twelve)
    ia = rep(-log(n), n)
    ip = - gamma * (r - r1) # check this
    z = ia + ip
    chiphi = rep(0, n)
    for(j in (1:n)){
      X = 0
      Y = 0
      rj= r[j]
      deltaj = delta[j]
      for(k in (1:n)[-j]){
        b = beta[k,j]
        rk = r[k]
        deltak = delta[k]
        denom = (deltaj + deltak) * (b + deltak)
        # lemma 1
        if(rj < rk){
          w = deltaj / denom  * exp(- deltak * (rk - (rj - lag)))
          x = deltak * w
          y = 1 - b * w
        } else{
          w = deltak / denom * exp(- deltaj * ((rj - lag) - rk))
          x = deltaj * w
          y = deltak / (b + deltak) + b * w
        }
        # line twelve
        X = X + b * x / y
        Y = Y + log(y)
      }
      chiphi[j] = Y + log(X)
    }
    # line eight
    for(alpha in 1:n){z[alpha] = z[alpha] + sum(chiphi[-alpha])}
    z = log(sum(exp(z)))
    a = sum(log(gamma / delta))
    # negative log likelihoods
    return(-(a+z))
  }
}
