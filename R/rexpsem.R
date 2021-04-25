#' Simulate Exponential Infectious Periods
#'
#' Draw exponential infectious periods assuming SEM model.
#'
#' @param beta numeric rate
#' @param gamma numeric rate
#' @param N integer population size
#'
#' @return infection and removal times
#'
#' @export
rexpsem = function(beta, gamma, N){
  # initialize vectors
  t = 0
  betaN = beta / N
  i = rep(Inf, N)
  r = rep(Inf, N)
  alpha = sample(N, 1)
  i[alpha] = t
  # simulate epidemic
  St = sum(is.infinite(i))
  It = sum(is.finite(i)) - sum(is.finite(r))
  while(It > 0){
    if(St > 0){
      # simulate infection times
      x = rexp(St, rate = betaN * It)
      minx = min(x)
      if(St > 1){
        argx = sample(which(is.infinite(i), arr.ind = T), 1)
      } else{
        argx = which(is.infinite(i), arr.ind = T)
      }
      # simulate removal times
      y = rexp(It, rate = gamma)
      miny = min(y)
      if(It > 1){
        argy = sample(which(is.infinite(r) & is.finite(i), arr.ind = T), 1)
      } else{
        argy = which(is.infinite(r) & is.finite(i))
      }
      # first arrival
      if(minx < miny){
        t = t + minx
        i[argx] = t
      } else{
        t = t + miny
        r[argy] = t
      }
    } else{
      # simulate removal times
      y = rexp(It, rate = gamma)
      miny = min(y)
      if(It > 1){
        argy = sample(which(is.infinite(r), arr.ind = T), 1)
      } else{
        argy = which(is.infinite(r), arr.ind = T)
      }
      t = t + miny
      r[argy] = t
    }
    # update (S,I) counts
    St = sum(is.infinite(i))
    It = sum(is.finite(i)) - sum(is.finite(r))
  }
  output = matrix(c(i,r),
                  nrow = N,
                  ncol = 2,
                  byrow = F)
  colnames(output) = c('i','r')
  return(output)
}
