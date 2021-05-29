#' Eichner-Dietz PBLA
#'
#' Compute the Eichner-Dietz likelihood approximation. Supports Erlang infectious periods.
#'
#' @param r numeric vector of increasing removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param m positive integer shape
#' @param nt integer points for trapezoidal integration
#' @param mint numeric lower bound for trapezoidal integration
#'
#' @return negative log likelihood
#'
#' @export
pbla_ed = function(r, beta, gamma, m = 1, nt = 100, mint = -2){

  # copy and paste from caTools
  trapz = function (x, y){
    idx = 2:length(x)
    return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
  }

  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }

  if((any(beta < 0)) | (any(gamma < 0)) | (!is.wholenumber(m)) | (m <= 0)){
    # invalid parameters
    return(1e15)
  } else{

    if(m == 1){
      # exponential infectious periods

      # initialize
      n = length(r)
      N = ncol(beta)
      t = seq(mint, max(r), length.out = nt)
      A = rep(0, nt)
      Y = 0

      # evaluate integrals
      for(j in 1:n){
        gammaj = gamma[j]
        rj = r[j]
        # calculate A
        for(k in 1:nt){
          A[k] = sum((beta[1:n,j] / gamma * exp(- gamma * (r - pmin(t[k], r))))[-j])
        }
        # calculate values for trapezium rule
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
    } else{
      # Erlang case

      # initialize
      n = length(r)
      N = ncol(beta)
      t = seq(mint, max(r), length.out = nt)

      #store some sum values
      U = matrix(0, nrow = n, ncol = nt)
      V = matrix(0, nrow = n, ncol = nt)
      for(j in 1:n){
        rj = r[j]
        gammaj = gamma[j]
        for(k in 1:nt){
          u = 0
          v = 0
          tk = t[k]
          for(l in 0:(m-1)){
            u = u + ((gammaj * (rj - min(tk, rj))) ^ l) / factorial(l) * (m - l)
            v = v + ((gammaj * (rj - tk)) ^ l) / factorial(l)
          }
          U[j,k] = u
          V[j,k] = v
        }
      }

      # evaluate integrals
      C = rep(0, nt)
      Y = 0
      for(j in 1:n){
        gammaj = gamma[j]
        rj = r[j]
        # calculate C
        for(k in 1:nt){
          C[k] = sum((beta[1:n,j] / gamma * exp(- gamma * (r - pmin(t[k], r))) * U[1:n,k])[-j])
        }
        # calculate values for trapezium rule
        y = rep(0, nt)
        for(k in (1:n)[-j]){
          rk = r[k]
          indices = which(t < min(rj, rk))
          y[indices] = y[indices] +
            beta[k,j] * exp(- gamma[k] * (rk - t[indices]) - C[indices]) * V[k,indices]
        }
        y = y * dgamma(rj - t, m, gammaj)
        Y = Y + log(trapz(t, y))
      }

      # failure to infect non-infectives
      Z = 0
      for(j in (n+1):N){Z = Z + sum(beta[1:n,j] * m / gamma)}

      # negative log likelihood
      return(Z - Y)
    }
  }
}
