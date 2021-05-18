# PBLA MCMC Example
# Seth Temple, sdtemple@uw.edu

library(pblas)

# Import Epidemic ---------------------------------------------------------

epidemic = readRDS("some-epidemic.rds")
r = epidemic[,2][is.finite(epidemic[,2])]
N = nrow(epidemic)
n = length(r)

# Set Parameter Values ----------------------------------------------------

K = 500 # iteration
B = 400 # burn-in
U = 100 # print frequency

pbla = pbla_std_gsem
m = 1 # Erlang shape

bg = nlm(pbla_gsem, c(1,1), pbla, r, N, c(m, n))$estimate

bsd = 0.2 # beta random walk standard deviation
gsd = 0.2 # gamma random walk standard deviation
binit = bg[1] # initial beta value
ginit = bg[2] # initial gamma value
bshape = 1 # beta prior shape
brate = 1e-4 # beta prior rate
gshape = 1 # gamma prior shape
grate = 1e-4 # gamma prior rate

# MCMC --------------------------------------------------------------------

#initialize
b = binit
g = ginit
storage = array(NA, dim = c(2, K))
ell = pbla(r, b, g, N, m)

# sampling
for(k in 1:K){
  
  # beta metropolis step
  bp = rnorm(1, b, bsd)
  if(bp > 0){
    ellp = - pbla(r, bp, g, N, m, n)
    a = min(1, exp(ellp + ell +
                     log(dgamma(bp, bshape, brate)) -
                     log(dgamma(b, bshape, brate))))
    if(runif(1) < a){
      b = bp
      ell = - ellp
    }
  }
  storage[1,k] = b
  
  # gamma metropolis step
  gp = rnorm(1, g, gsd)
  if(gp > 0){
    ellp = - pbla(r, b, gp, N, m, n)
    a = min(1, exp(ellp + ell +
                     log(dgamma(gp, gshape, grate)) -
                     log(dgamma(g, gshape, grate))))
    if(runif(1) < a){
      g = gp
      ell = - ellp
    }
  }
  storage[2,k] = g
  
  # print iteration
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
