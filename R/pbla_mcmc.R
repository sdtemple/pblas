#' PBLA MCMC for Stochastic Epidemic Models
#'
#' Use PBLA in MCMC with Metropolis-Hastings-based random walk.
#'
#'  @param niter integer iterations
#'  @param sd standard deviation of random walk proposal
#'  @param pbla likelihood approximation
#'  @param betamap matrix of assignments
#'  @param gammamap vector of assignments
#'  @param r numeric vector of increasing removal times
#'  @param etc other parameters to pass (e.g. lag)
#'  @param betastart positive initial value
#'  @param betashape positive prior
#'  @param betarate positive prior
#'  @param gammastart positive initial value
#'  @param gammashape positive prior
#'  @param gammarate positive prior
#'
#'  @return posterior sample
#'
#'  @export
pbla_mcmc = function(niter, sd,
                     pbla, betamap, gammamap, r, etc,
                     betastart, betashape, betarate,
                     gammastart, gammashape, gammarate){

  storage = array(NA, dim = c(2, niter))
  n = length(r)
  beta = betastart
  gamma = gammastart
  ell = pbla_het(c(beta, gamma), 1, pbla, betamap, gammamap, r, etc)
  for(k in 1:niter){

    # beta metropolis step
    betap = rnorm(1, beta, sd)
    if(betap > 0){
      ellp = - pbla_het(c(betap, gamma), 1, pbla, betamap, gammamap, r, etc)
      a = min(1, exp(ellp + ell -
                       log(dgamma(betap, betashape, betarate)) +
                       log(dgamma(beta, betashape, betarate))))
      if(runif(1) < a){
        beta = betap
        ell = - ellp
      }
    }
    storage[1,k] = beta

    # gamma metropolis step
    gammap = rnorm(1, gamma, sd)
    if(gammap > 0){
      ellp = - pbla_het(c(beta, gammap), 1, pbla, betamap, gammamap, r, etc)
      a = min(1, exp(ellp + ell -
                       log(dgamma(gammap, gammashape, gammarate)) +
                       log(dgamma(gamma, gammashape, gammarate))))
      if(runif(1) < a){
        gamma = gammap
        ell = - ellp
      }
    }
    storage[2,k] = gamma

    if(!(k %% 1000)){
      print(k)
    }

  }

  return(storage)
}
