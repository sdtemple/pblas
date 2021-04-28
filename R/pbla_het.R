#' PBLAs for Stochastic Epidemic Model (Heterogeneous Mixing)
#'
#' Run a pair-based likelihood approximation for a stochastic epidemic model. Heterogeneous mixing available via mapping assignments.
#'
#' @param rates numeric vector of rates
#' @param R integer index to split rates on
#' @param pbla function
#' @param betamap matrix of assignments
#' @param gammamap vector of assignments
#' @param r numeric vector of increasing removal times
#' @param etc other parameters to pass (e.g. lag)
#'
#' @return negative log likelihood
#'
#' @export
pbla_het = function(rates, R, pbla, betamap, gammamap, r, etc){
  br = rates[1:R]
  gr = rates[(R+1):length(rates)]
  beta = map_beta(br, betamap)
  beta = beta / ncol(beta)
  gamma = gr[gammamap]
  return(do.call(pbla, c(list(r=r,beta=beta,gamma=gamma), etc)))
}
