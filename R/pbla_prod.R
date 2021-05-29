#' Product PBLA
#'
#' Based on product independence, compute pair-based likelihood approximation. Assumes exponential infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta numeric rate
#' @param gamma numeric rate
#' @param N integer population size
#' @param A integer patient zeros
#'
#' @return negative log likelihood
#'
#' @export
pbla_prod = function(r, beta, gamma, N, A = 1){

  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }

  if((any(beta < 0)) | (any(gamma < 0)) |
     (!is.wholenumber(N)) | (N <= 0) |
     (!is.wholenumber(A)) | (A <= 0)){
    # invalid parameters
    return(1e15)
  } else{

    # initialize
    n = length(r)
    r1 = r[1]
    beta = beta / N
    lb = log(beta)
    l2 = log(2)

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
    pe = sum(log(delta) - log(beta * 1:(n-1) + delta)) # product expectation (lemma 2)
    z = ia + ip + pe

    # evaluate chi terms
    for(j in (1:n)){
      X = 0
      rj = r[j]
      for(k in (1:n)[-j]){
        rk = r[k]
        if(rj < rk){
          X = X + exp(- delta * (rk - rj))
        } else{
          X = X + exp(- delta * (rj - rk))
        }
      }
      z[-j] = z[-j] + log(X) + lb - l2
    }

    # line 8
    z = matrixStats::logSumExp(z)
    a = n * (log(gamma) - log(delta))

    # negative log likelihood
    return(-(a+z))
  }
}
