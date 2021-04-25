#' PBLAs for Stochastic Epidemic Model (Homogeneous Mixing)
#'
#' Run a pair-based likelihood approximation for a stochastic epidemic model. Heterogeneous mixing available via mapping assignments.
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
pbla_homo = function(rates, pbla, r, N, etc){
  beta = rates[1]
  gamma = rates[2]
  return(do.call(pbla, c(list(r=r,beta=beta,gamma=gamma,N=N), etc)))
}
