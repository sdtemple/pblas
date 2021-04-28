#' Non-Centered Data-Augmented MCMC Sampler for SEM
#'
#' Draw posterior samples for a stochastic epidemic model.
#'
#' @param niter integer iterations
#' @param sd positive control of random walk proposal
#' @param J integer infection time updates
#' @param r numeric vector of removal times
#' @param N integer population size
#' @param m integer Erlang shape
#' @param betastart positive
#' @param betashape prior
#' @param betarate prior
#' @param gammastart positive
#' @param gammashape prior
#' @param gammarate prior
#'
#' @return array
#'
#' @export
ncda_mcmc_sem = function(niter, sd, J, # mcmc parameters 
                         r, N, m, # real data
                         betastart, betashape, betarate, # beta priors 
                         gammastart, gammashape, gammarate){ # gamma priors
  
  # partial log likelihood
  pell_pi_gamma_i = function(r, i, gamma, m, betashape, betarate, gammashape, gammarate){
    
    # initialize
    n = length(r)
    N = length(i)
    tau = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      tau[j,] = sapply(i, min, r[j]) - sapply(i, min, i[j])
    }
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < i[j]) * (r > i[j])
    }
    Blong = apply(ind, 1, sum)
    Blong = Blong[Blong > 0]
    
    # log likelihood
    ell = log(dgamma(gamma, shape = gammashape, rate = gammarate))
    ell = ell + sum(log(Blong))
    ell = ell + sum(log(dgamma(r - i[1:n], shape = m, rate = gamma)))
    ell = ell - (betashape + n - 1) * log(betarate + sum(tau))
    return(ell)
  }
  
  ell_hastings_i = function(r, i, ip, betashape, betarate){
    # initialize
    n = length(r)
    N = length(i)
    tau = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      tau[j,] = sapply(i, min, r[j]) - sapply(i, min, i[j])
    }
    taup = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      taup[j,] = sapply(ip, min, r[j]) - sapply(ip, min, ip[j])
    }
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < i[j]) * (r > i[j])
    }
    Blong = apply(ind, 1, sum)
    Blong = Blong[Blong > 0]
    indp = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      indp[j,] = (ip[1:n] < ip[j]) * (r > ip[j])
    }
    Bplong = apply(indp, 1, sum)
    Bplong = Bplong[Bplong > 0]
    
    ellratio = sum(log(Bplong)) - sum(log(Blong))
    ellratio = ellratio + (betashape + n - 1) * 
      (log(betarate + sum(tau)) - log(betarate + sum(taup)))
    return(ellratio)
  }
  
  # initialize
  storage = array(NA, dim = c(2, niter))
  n = length(r)
  tau = matrix(0, nrow = n, ncol = N)
  #taup = matrix(0, nrow = n, ncol = N)
  
  beta = betastart / N
  gamma = gammastart
  
  i = r - (rgamma(n, shape = m, rate = 1) / gamma)
  i = c(i, rep(Inf, N - n))
  while(!is_epidemic(r, i)){ # must be epidemic
    i = r - (rgamma(n, shape = m, rate = 1) / gamma)
    i = c(i, rep(Inf, N - n))
  }
  
  # mcmc sampling
  for(k in 1:niter){
    
    # gamma metropolis hastings step
    gammap = rnorm(1, mean = gamma, sd = sd)
    if(gammap > 0){
      ip = r - ((r - i[1:n]) * gamma / gammap)
      ip = c(ip, rep(Inf, N - n))
      if(is_epidemic(r, ip)){
        a = min(1, exp(pell_pi_gamma_i(r, ip, gammap, m, betashape, betarate, gammashape, gammarate) -
                         pell_pi_gamma_i(r, i, gamma, m, betashape, betarate, gammashape, gammarate)))
        if(runif(1) < a){
          gamma = gammap
          i = ip
        }
      }
    #  for(j in 1:n){
    #    taup[j,] = sapply(ip, min, r[j]) - sapply(ip, min, ip[j])
    #  }
    #  betap = rgamma(1, shape = betashape + n - 1, rate = betarate + sum(taup))
    #  a = min(1, exp(aloglike_sem(r, ip, betap, gammap, m) - 
    #            aloglike_sem(r, i, beta, gamma, m)))
    #  if((runif(1) < a)){
    #    beta = betap
    #    gamma = gammap
    #    i = ip
    #  }
    }
    storage[2,k] = gamma
    
    # infection times metropolis hastings step
    for(j in 1:J){
      l = sample(1:n, 1)
      il = r[l] - (rgamma(1, shape = m, rate = 1) / gamma)
      ip = i
      ip[l] = il
      if(is_epidemic(r, ip)){ # must be epidemic
        a = min(1, exp(ell_hastings_i(r, i, ip, betashape, betarate)))
        #a = min(1, exp(aloglike_sem(r, ip, beta, gamma, m) -
        #                 aloglike_sem(r, i, beta, gamma, m) +
        #                 log(dgamma(r[j] - i[j], shape = m,  rate = gamma)) -
        #                 log(dgamma(r[j] - ip[j], shape = m, rate = gamma))))
        if(runif(1) < a){
          i[l] = il
        }
      }
    }
    
    # beta gibbs step
    for(j in 1:n){
      tau[j,] = sapply(i, min, r[j]) - sapply(i, min, i[j])
    }
    beta = rgamma(1, shape = betashape + n - 1, rate = betarate + sum(tau))
    storage[1,k] = beta * N
    
    # iteration update
    if(!(k %% 100)){
      print(k)
    }
  }
  
  return(storage)
}
