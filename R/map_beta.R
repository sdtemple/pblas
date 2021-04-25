#' Assign Infection Rates
#'
#' Assign individual infection rates.
#'
#' @param beta numeric vector of infection rates
#' @param map matrix of assignments
#'
#' @return matrix of infection rates
#'
#' @export
map_beta = function(beta, map){
  n = nrow(map)
  N = ncol(map)
  betamap = matrix(beta[map], nrow = n, ncol = N)
  return(betamap)
}
