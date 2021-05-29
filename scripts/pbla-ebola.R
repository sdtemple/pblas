# PBLA Ebola Virus in West Africa
# Seth Temple, sdtemple@uw.edu
# May 17, 2021

# Data Prep ---------------------------------------------------------------

guinea = readRDS("guinea-deaths.rds")
sierra = readRDS("sierra-deaths.rds")
liberia = readRDS("liberia-deaths.rds")

# initialize
scalar = 1e3
gamma = 1 / 5.61 * scalar # removal rate
lag = 5.3 / scalar # exposed to infectious rate
N = 1e6 # population size
A = 1 # patient zeros

# Custom Functions --------------------------------------------------------

# simplified version of pbla_std
pbla = function(r){

  if(any(beta < 0)){
    return(1e15)
  }

  r1 = r[1]
  delta = B + gamma

  # calculate log likelihood (line six)
  ia = rep(-log(A), A)
  ip = - delta[1:A] * (r[1:A] - r1)
  z = ia + ip

  # evaluate psi and chi terms
  psichi = rep(0, n)
  for(j in (1:n)){
    X = 0
    Y = 0
    rj= r[j]
    deltaj = delta[j]
    for(k in (1:n)[-j]){
      b = beta[k,j]
      rk = r[k]
      deltak = delta[k]
      denom = (deltaj + deltak) * (b + deltak)
      # lemma 1
      if(rj - lag < rk){
        w = deltaj / denom  * exp(- deltak * (rk - (rj - lag)))
        x = deltak * w
        y = 1 - b * w
      } else{
        w = deltak / denom * exp(- deltaj * ((rj - lag) - rk))
        x = deltaj * w
        y = deltak / (b + deltak) + b * w
      }
      # line twelve
      X = X + b * x / y
      Y = Y + log(y)
    }
    psichi[j] = Y + log(X)
  }

  # line eight
  for(alpha in 1:A){z[alpha] = z[alpha] + sum(psichi[-alpha])}
  z = matrixStats::logSumExp(z)
  a = sum(log(gamma) - log(delta))

  # negative log likelihoods
  return(-(a+z))
}

# use <<- to update large global beta matrix
# use <<- to update large global B vector
pbla_ebola_seir = function(x, r){
  b0 = x[1]
  k0 = x[2]
  n = length(r)
  for(j in 1:n){
    rj = r[j]
    for(k in (1:n)[-j]){
      rk  = r[k]
      if(rj > rk - lag){
        tjk = rk - (1 / gamma) - lag - (exp(- gamma * (rj - rk + lag)) / 4 / gamma)
        beta[j,k] <<- b0 * exp(- k0 * tjk) / N #global update
      } else{
        tjk = rj - (1 / 2 / gamma) - (3 * exp(- gamma * (rk - rj - lag)) / 4 / gamma)
        beta[j,k] <<- b0 * exp(- k0 * tjk) / N # global update
      }
    }
  }
  for(j in 1:n){
    tjk = r[j] - (1 / 2 / gamma)
    B[j] <<- (N - n) * b0 * exp(- k0 * tjk) / N # global update
  }
  return(pbla(r))
}

# Guinea ----------------------------------------------------------------

r = guinea / scalar
n = length(r)
beta = matrix(0, n, n)
B = rep(0, n)
pbla_ebola_seir(c(.2,.001)*scalar, r)

b0init = 0.4
k0init = 0.002
ptm = proc.time()
out = nlm(pbla_ebola_seir, c(b0init, k0init)*scalar, r,
          print.level = 2,
          hessian = T)
saveRDS(out, "guinea-nlm.rds")
ptm = (proc.time() - ptm)[3]
print(ptm / 60)
# save these

# Sierra Leone ------------------------------------------------------------

r = sierra / scalar
n = length(r)
beta = matrix(0, n, n)
B = rep(0, n)

b0init = 1
k0init = 0.01
ptm = proc.time()
out = nlm(pbla_ebola_seir, c(b0init, k0init)*scalar, r,
          print.level = 2,
          hessian = T)
saveRDS(out, "sierra-nlm.rds")
ptm = (proc.time() - ptm)[3]
print(ptm / 60)
# save these

# Liberia ------------------------------------------------------------

r = liberia / scalar
n = length(r)
beta = matrix(0, n, n)
B = rep(0, n)

b0init = 2
k0init = 0.05
ptm = proc.time()
out = nlm(pbla_ebola_seir, c(b0init, k0init)*scalar, r,
          print.level = 2,
          hessian = T)
saveRDS(out, "liberia-nlm.rds")
ptm = (proc.time() - ptm)[3]
print(ptm / 60)
# save these
