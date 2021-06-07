#' MLE for Completely Observed SEM
#'
#' Compute MLE for completely observed stochastic epidemic.
#'
#' @param r numeric vector of removal times
#' @param i numeric vector of infection times
#' @param N integer population size
#'
#' @return MLE for (beta, gamma)
#'
#' @export
mle_gsem = function(r, i, N){
    n = length(r)
    t = 0
    for(j in 1:n){
        t = t + sum(sapply(i, min, r[j]) - sapply(i, min, i[j]))
    }
    ri = sum(r - i)
    g = ri / n
    b = (n - 1) / (t + (N - n) * ri) * N
    return(c(b,g))
}