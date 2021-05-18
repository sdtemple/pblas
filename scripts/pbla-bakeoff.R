# PBLA Bake Off
# Seth Temple, sdtemple@uw.edu

library(pblas)

# Set Parameter Values ----------------------------------------------------------

# epidemic
N = 100
A = 5
beta = 1.5
gamma = 1

# file name
fn = paste("b", beta, "-g", gamma, "-N", N,
           "-bakeoff", index, ".rds")

# Simulation Study ----------------------------------------------------------

# initialize
ct = 0
K = 100
U = 10
storage = array(NA, dim = c(11,K))

# run
while(ct < K){
  
  epi = rgsem(beta, gamma, N)
  r = (epi[,2])[is.finite(epi[,2])]
  r = sort(r)
  n = length(r)
  while(n < A){
    epi = rgsem(beta, gamma, N)
    r = (epi[,2])[is.finite(epi[,2])]
    r = sort(r)
    n = length(r)
  }
  
  ct = ct + 1
  storage[1,ct] = n
  storage[2:3,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_std_gsem, r=r, N=N, c(1, n))$estimate
  storage[4:5,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_prod, r=r, N=N, n)$estimate
  storage[6:7,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_weak, r=r, N=N, n)$estimate
  storage[8:9,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_sep_gsem, r=r, N=N, n)$estimate
  storage[10:11,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_f_gsem, r=r, N=N, n)$estimate
  
  if(!(ct %% U)){
    print(ct)
  }
  
}
saveRDS(fn, storage)