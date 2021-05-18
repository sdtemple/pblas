# PBLA Simulation Study (Erlang)
# Seth Temple, sdtemple@uw.edu

library(pblas)

# Set Parameter Values ----------------------------------------------------------

beta = 1.5
gamma = 1
N = 100
m = 2

# epidemic
A = 5

# file name
fn = paste("b", beta, "-g", gamma, "-N", N, "-m", m, 
           ".rds", sep = "")

# Simulation Study ----------------------------------------------------------

# initialize
ct = 0
K = 50
U = 10
storage = array(NA, dim = c(5,K))

# run
while(ct < K){
  
  epi = rgsem(beta, gamma, N, m)
  r = (epi[,2])[is.finite(epi[,2])]
  r = sort(r)
  n = length(r)
  while(n < A){
    epi = rgsem(beta, gamma, N, m)
    r = (epi[,2])[is.finite(epi[,2])]
    r = sort(r)
    n = length(r)
  }
  
  ct = ct + 1
  storage[1,ct] = n
  storage[2:3,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_std_gsem, r=r, N=N, m)$estimate
  storage[4:5,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_ed_gsem, r=r, N=N, m)$estimate
  
  if(!(ct %% U)){
    print(ct)
  }
  
}
saveRDS(storage, fn)