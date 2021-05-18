# PBLA MCMC for Tristan da Cunha
# Seth Temple, sdtemple@uw.edu

library(pblas)

# Import Epidemic ---------------------------------------------------------

cases = read.csv("tdc_jitteredtimes.txt", header = F)
ages = read.csv("tdc_agegroups.txt", header = F)
tdc = cbind(cases, ages)
colnames(tdc) = c('day','age')
tdc$day[tdc$day == -1000] = Inf

N = nrow(tdc) # population size
N1 = sum(tdc$age == 1)
N2 = sum(tdc$age == 2)
N3 = sum(tdc$age == 3)
n = sum(is.finite(tdc$day)) # infecteds
r = tdc$day[1:n]

# beta map
bm = matrix(0, nrow = N, ncol = N)
for(j in 1:N){
  bm[,j] = tdc$age[j]
}

# gamma map
gm = rep(1, n)


# Set Parameter Values ----------------------------------------------------

set.seed(5152021)

K = 11000 # iteration
U = 1000 # print frequency

bsd = 0.005 # beta random walk standard deviation
gsd = 0.2 # gamma random walk standard deviation

b1init = 0.001
b2init = 0.001
b3init = 0.001
ginit = 0.5

bshape = 1e-8 # beta prior shape
brate = 1e-5 # beta prior rate
gshape = 1e-4 # gamma prior shape
grate = 1e-3 # gamma prior rate

pbla = pbla_std # standard pair-based likelihood approximation

# MCMC --------------------------------------------------------------------

#initialize
b1 = b1init
b2 = b2init
b3 = b3init
g = ginit
storage = array(NA, dim = c(5, K))
ell = pbla_multi(c(b1, b2, b3, g), 3, pbla, bm, gm, r)

for(k in 1:K){

  # b1 metropolis step
  b1p = rnorm(1, b1, bsd)
  if(b1p > 0){
    ellp = - pbla_multi(c(N * c(b1p, b2, b3), g), 3, pbla, bm, gm, r)
    a = min(1, exp(ellp + ell +
                     log(dgamma(b1p, bshape, brate)) -
                     log(dgamma(b1, bshape, brate))))
    if(runif(1) < a){
      b1 = b1p
      ell = - ellp
    }
  }
  storage[1,k] = b1

  # b2 metropolis step
  b2p = rnorm(1, b2, bsd)
  if(b2p > 0){
    ellp = - pbla_multi(c(N * c(b1, b2p, b3), g), 3, pbla, bm, gm, r)
    a = min(1, exp(ellp + ell +
                     log(dgamma(b2p, bshape, brate)) -
                     log(dgamma(b2, bshape, brate))))
    if(runif(1) < a){
      b2 = b2p
      ell = - ellp
    }
  }
  storage[2,k] = b2

  # b3 metropolis step
  b3p = rnorm(1, b3, bsd)
  if(b3p > 0){
    ellp = - pbla_multi(c(N * c(b1, b2, b3p), g), 3, pbla, bm, gm, r)
    a = min(1, exp(ellp + ell +
                     log(dgamma(b3p, bshape, brate)) -
                     log(dgamma(b3, bshape, brate))))
    if(runif(1) < a){
      b3 = b3p
      ell = - ellp
    }
  }
  storage[3,k] = b3

  # g metropolis step
  gp = rnorm(1, g, gsd)
  if(gp > 0){
    ellp = - pbla_multi(c(N * c(b1, b2, b3), gp), 3, pbla, bm, gm, r)
    a = min(1, exp(ellp + ell +
                     log(dgamma(gp, gshape, grate)) -
                     log(dgamma(g, gshape, grate))))
    if(runif(1) < a){
      g = gp
      ell = - ellp
    }
  }
  storage[4,k] = g

  storage[5,k] = (storage[1,k] * N1 + storage[2,k] * N2 + storage[3,k] * N3) / storage[4,k]

  # print iteration
  if(!(k %% U)){
    print(k)
  }

}

# saving
saveRDS(storage, "pbla-mcmc-tdc.rds")


# Analysis ------------------------------------------------------

# traceplots
par(mfrow=c(2,2))
plot(storage[1,], type = "l", ylab = expression(beta[1]))
plot(storage[2,], type = "l", ylab = expression(beta[2]))
plot(storage[3,], type = "l", ylab = expression(beta[3]))
plot(storage[4,], type = "l", ylab = expression(gamma))

# posterior means
burnin = 1000
pb1 = mean(storage[1,burnin:K])
pb2 = mean(storage[2,burnin:K])
pb3 = mean(storage[3,burnin:K])
pg = mean(storage[4,burnin:K])
pr = mean(storage[5,burnin:K])

# Hayakawa et al. (2003)
hb1 = 0.00451
hb2 = 0.00181
hb3 = 0.00131
hg = 0.371
hr = (hb1 * N1 + hb2 * N2 + hb3 * N3) / hg

# pbla mle
mle = nlm(pbla_multi, c(N * c(b1init, b2init, b3init), ginit),
          R=3, pbla=pbla_std, betamap=bm, gammamap=gm, r=r)$estimate
mle = c(mle[1:3]/N, mle[4])
mle = c(mle, (mle[1] * N1 + mle[2] * N2 + mle[3] * N3) / mle[4])

# ed mle
ed = nlm(pbla_multi, c(N * c(b1init, b2init, b3init), ginit),
         R=3, pbla=pbla_ed, betamap=bm, gammamap=gm, r=r,
         c(1, 1000, -5))$estimate # wait a couple minutes
ed = c(ed[1:3]/N, ed[4])
ed = c(ed, (ed[1] * N1 + ed[2] * N2 + ed[3] * N3) / ed[4])

# Table -------------------------------------------------------------------

tbl = matrix(NA, 5, 4)
colnames(tbl) = c('PBLA MCMC', 'DAMCMC', 'E+D MLE', 'PBLA MLE')

# pbla mcmc
tbl[1,1] = round(pb1,5)
tbl[2,1] = round(pb2,5)
tbl[3,1] = round(pb3,5)
tbl[4,1] = round(pg,5)
tbl[5,1] = round(pr,5)

# damcmc
tbl[1,2] = round(hb1,5)
tbl[2,2] = round(hb2,5)
tbl[3,2] = round(hb3,5)
tbl[4,2] = round(hg,5)
tbl[5,2] = round(hr,5)

# e+d mle
tbl[1,3] = round(ed[1],5)
tbl[2,3] = round(ed[2],5)
tbl[3,3] = round(ed[3],5)
tbl[4,3] = round(ed[4],5)
tbl[5,3] = round(ed[5],5)

# pbla mle
tbl[1,4] = round(mle[1],5)
tbl[2,4] = round(mle[2],5)
tbl[3,4] = round(mle[3],5)
tbl[4,4] = round(mle[4],5)
tbl[5,4] = round(mle[5],5)

stargazer::stargazer(tbl, digits = 5)

# Plots -------------------------------------------------------------------

acol = c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
acol = c(acol[6], acol[7], acol[4])

par(mfrow=c(2,2))
hist(storage[1,burnin:K],
     breaks = 30,
     main = NA,
     xlab = expression(beta[1]),
     xlim = c(0,0.015),
     col = "white")
abline(v = hb1, col = acol[1], lty = 2, lwd = 2)
abline(v = ed[1], col = acol[2], lty = 3, lwd = 2)
abline(v = mle[1], col = acol[3], lty = 4, lwd = 2)
legend("topright",
       c("DAMCMC","E+D MLE", "PBLA MLE"),
       col = acol,
       lwd = rep(2, 3),
       lty = 2:4)

hist(storage[2,burnin:K],
     breaks = 30,
     main = NA,
     xlab = expression(beta[2]),
     xlim = c(0,0.006),
     col = "white")
abline(v = hb2, col = acol[1], lty = 2, lwd = 2)
abline(v = ed[2], col = acol[2], lty = 3, lwd = 2)
abline(v = mle[2], col = acol[3], lty = 4, lwd = 2)

hist(storage[3,burnin:K],
     breaks = 30,
     main = NA,
     xlab = expression(beta[3]),
     xlim = c(0,0.004),
     col = "white")
abline(v = hb3, col = acol[1], lty = 2, lwd = 2)
abline(v = ed[3], col = acol[2], lty = 3, lwd = 2)
abline(v = mle[3], col = acol[3], lty = 4, lwd = 2)

hist(storage[4,burnin:K],
     breaks = 30,
     main = NA,
     xlab = expression(gamma),
     xlim = c(0,1),
     col = "white")
abline(v = hg, col = acol[1], lty = 2, lwd = 2)
abline(v = ed[4], col = acol[2], lty = 3, lwd = 2)
abline(v = mle[4], col = acol[3], lty = 4, lwd = 2)
