#' PBLAs for Stochastic Epidemic Model (Heterogeneous Mixing)
#'
#' Run a pair-based likelihood approximation for a stochastic epidemic model. Heterogeneous mixing available via mapping assignments. Compatible with `pbla_std`, `pbla_ed`, `pbla_sep`, and `pbla_f`.
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
pbla_het = function(rates, R, pbla, betamap, gammamap, r, etc = NULL){
  br = rates[1:R]
  gr = rates[(R+1):length(rates)]
  beta = map_beta(br, betamap)
  beta = beta / ncol(beta)
  gamma = gr[gammamap]
  if(is.null(etc)){ # use defaults
    return(do.call(pbla, list(r=r,beta=beta,gamma=gamma)))
  } else{ # pass in all parameters
    return(do.call(pbla, c(list(r=r,beta=beta,gamma=gamma), etc)))
  }
}
