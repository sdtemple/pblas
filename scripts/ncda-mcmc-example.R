# Non-Centered DAMCMC Example
# Seth Temple, sdtemple@uw.edu

library(pblas)

# Import Epidemic ---------------------------------------------------------

epidemic = readRDS("some-epidemic.rds")
r = epidemic[,2][is.finite(epidemic[,2])]
n = length(r)
N = nrow(epidemic)

# Set Parameter Values ----------------------------------------------------

K = 500 # iteration
B = 400 # burn-in
J = 2 # per iteration infection time updates
U = 100 # print frequency

m = 1 # Erlang shape

bg = nlm(pbla_gsem, c(1,1), pbla_std_gsem, r, N, c(m, n))$estimate

gsd = 0.2 # gamma random walk standard deviation
binit = bg[1] # initial beta value
ginit = bg[2] # initial gamma value
bshape = 1 # beta prior shape
brate = 1e-4 # beta prior rate
gshape = 1 # gamma prior shape
grate = 1e-4 # gamma prior rate

# Custom Functions --------------------------------------------------------

# indicates data consistent with epidemic
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

# utility for recovery rate metropolis hastings step
gupdate = function(r, i, g, m, bshape, brate, gshape, grate){
  
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
  ell = log(dgamma(g, shape = gshape, rate = grate))
  ell = ell + sum(log(Blong))
  ell = ell - (bshape + n - 1) * log(brate + sum(tau))
  return(ell)
}

# utlity for infection time metropolis hastings step
iupdate = function(r, i, ip, bshape, brate){
  
  # initialize
  n = length(r)
  N = length(i)
  
  # compute tau matrices
  tau = matrix(0, nrow = n, ncol = N)
  for(j in 1:n){
    tau[j,] = sapply(i, min, r[j]) - sapply(i, min, i[j])
  }
  taup = matrix(0, nrow = n, ncol = N)
  for(j in 1:n){
    taup[j,] = sapply(ip, min, r[j]) - sapply(ip, min, ip[j])
  }
  
  # compute 
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
  ellratio = ellratio + (bshape + n - 1) *
    (log(brate + sum(tau)) - log(brate + sum(taup)))
  return(ellratio)
}

# MCMC --------------------------------------------------------------------

# initialize
storage = array(NA, dim = c(2, K))
tau = matrix(0, nrow = n, ncol = N)

b = binit / N
g = ginit

i = r - (rgamma(n, shape = m, rate = 1) / g)
i = c(i, rep(Inf, N - n))
while(!is_epidemic(r, i)){ # must be epidemic
  i = r - (rgamma(n, shape = m, rate = 1) / g)
  i = c(i, rep(Inf, N - n))
}

# sampling
for(k in 1:K){
  
  # gamma metropolis hastings step
  gp = rnorm(1, mean = g, sd = gsd)
  if(gp > 0){
    ip = r - ((r - i[1:n]) * g / gp)
    ip = c(ip, rep(Inf, N - n))
    if(is_epidemic(r, ip)){
      a = min(1, exp(gupdate(r, ip, gp, m, bshape, brate, gshape, grate) -
                       gupdate(r, i, g, m, bshape, brate, gshape, grate)))
      if(runif(1) < a){
        g = gp
        i = ip
      }
    }
  }
  storage[2,k] = g
  
  # infection times metropolis hastings step
  for(j in 1:J){
    l = sample(1:n, 1)
    il = r[l] - (rgamma(1, shape = m, rate = 1) / g)
    ip = i
    ip[l] = il
    if(is_epidemic(r, ip)){ # must be epidemic
      a = min(1, iupdate(r, i, ip, bshape, brate))
      if(runif(1) < a){
        i[l] = il
      }
    }
  }
  
  # beta gibbs step
  for(j in 1:n){
    tau[j,] = sapply(i, min, r[j]) - sapply(i, min, i[j])
  }
  b = rgamma(1, shape = bshape + n - 1, rate = brate + sum(tau))
  storage[1,k] = b * N
  
  # iteration update
  if(!(k %% U)){
    print(k)
  }
  
  if(k == B){
    ptm = proc.time()
  }
}
t = unname((proc.time() - ptm)[3])

# Trace Plots -------------------------------------------------------------

par(mfrow=c(1,2))
plot(storage[1,], type = 'l')
plot(storage[2,], type = 'l')
print(mean(storage[1,B:K]))
print(mean(storage[2,B:K]))

# Effective Sample Size ---------------------------------------------------

be = LaplacesDemon::ESS(storage[1,B:K])
ge = LaplacesDemon::ESS(storage[2,B:K])
bet = log(be) - log(t)
get = log(ge) - log(t)
print(bet)
print(get)
