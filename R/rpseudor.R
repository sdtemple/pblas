#' Pseudo Removal Times
#'
#' Simulate new removal times from normal kernel density estimate.
#'
#' @param r removal times
#' @param eta proportion reported
#'
#' @return removal times
#'
#' @export
rpseudor = function(r, eta){
  n = length(r)
  nd = ceiling(n / eta) - n
  bw = density(r)$bw
  rn = rep(0, nd)
  for(j in 1:nd){
   rp = sample(r, 1)
   rn[j] = rnorm(1, rp, bw)
  }
  return(sort(c(r, rn)))
}
