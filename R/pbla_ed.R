#' Eichner-Dietz PBLA (Exponential)
#'
#' Compute the Eichner-Dietz likelihood approximation. Assume exponential infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param nt integer points for trapezoidal integration
#' @param mint numeric lower bound for trapezoidal integration
#'
#' @return negative log likelihood
#'
#' @export
pbla_ed = function(r, beta, gamma, nt, mint){
  if((any(beta <= 0)) | (any(gamma <= 0))){
    return(1e9) # positive rates
  } else{
    # copy and paste from caTools
    trapz = function (x, y){
      idx = 2:length(x)
      return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
    }
    # initialize
    n = length(gamma)
    t = seq(mint, max(r), length.out = nt)
    A = rep(0, nt)
    Y = 0
    # evaluate integrals
    for(j in 1:n){
      gammaj = gamma[j]
      rj = r[j]
      for(k in 1:nt){
        A[k] = sum((beta[1:n,j] / gamma * exp(- gamma * (r - pmin(t[k], r))))[-j])
      }
      y = rep(0, nt)
      for(k in (1:n)[-j]){
        rk = r[k]
        indices = which(t < min(rj, rk))
        y[indices] = y[indices] + 
          beta[k,j] * exp(- gamma[k] * (rk - t[indices]) - gammaj * (rj - t[indices]) - A[indices])  
      }
      Y = Y + log(gammaj) + log(trapz(t, y)) 
    }
    
    # failure to infect non-infectives
    Z = 0
    for(j in (n+1):N){Z = Z + sum(beta[1:n,j] / gamma)}
    # negative log likelihood
    return(Z - Y) 
  }
}
