#' PBLAs for Stochastic Epidemic Model (Homogeneous Mixing)
#'
#' Run a pair-based likelihood approximation for a stochastic epidemic model. Assumes homogeneous mixing. Compatible with `pbla_weak` and `pbla_prod`.
#'
#' @param rates numeric vector of rates
#' @param pbla function
#' @param r numeric vector of increasing removal times
#' @param N integer population size
#' @param etc other parameters to pass (e.g. lag)
#'
#' @return negative log likelihood
#'
#' @export
pbla_homo = function(rates, pbla, r, N, etc = NULL){
  beta = rates[1]
  gamma = rates[2]
  if(is.null(etc)){ # use defaults
    return(do.call(pbla, list(r=r,beta=beta,gamma=gamma,N=N)))
  } else{ # pass in all parameters
    return(do.call(pbla, c(list(r=r,beta=beta,gamma=gamma,N=N), etc)))
  }
}
