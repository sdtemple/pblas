#' Epidemic Comptability
#' 
#' Indicate if infection and removal times are compatible.
#' 
#' @param r numeric vector of removal times
#' @param i numeric vector of infection times
#' 
#' @return boolean
#' 
#' @export
is_epidemic = function(r, i){
  n = length(r)
  ind = matrix(0, nrow = n, ncol = n)
  for(j in 1:n){
    ind[j,] = (i[1:n] < i[j]) * (r > i[j])
  }
  x = apply(ind, 1, sum)
  x = x[x == 0]
  if(length(x) > 1){
    return(0)
  } else{
    return(1)
  }
}
