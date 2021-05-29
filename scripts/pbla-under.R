# PBLA Underreporting
# Seth Temple, sdtemple@uw.edu
# May 28, 2021

library(pblas)
set.seed(5212021)

N = 500
A = 30
beta = 1.5
gamma = 1

ct = 0
K = 500
U = 20
x = array(dim=c(11,K))

# run
while(ct < K){
  epi = rgsem(beta, gamma, N)
  r = (epi[,2])[is.finite(epi[,2])]
  r = sort(r)
  n = length(r)
  u = runif(1)
  keep = as.logical(rbinom(n, 1, u))
  keep = which(keep, arr.ind = T)
  ru = r[keep]
  nu = length(ru)
  while(nu < A){
    epi = rgsem(beta, gamma, N)
    r = (epi[,2])[is.finite(epi[,2])]
    r = sort(r)
    n = length(r)
    u = runif(1)
    keep = as.logical(rbinom(n, 1, u))
    keep = which(keep, arr.ind = T)
    ru = r[keep]
    nu = length(ru)
  }
  ct = ct + 1
  x[1,ct] = n
  x[2:3,ct] = nlm(pbla_gsem, c(2,1), pbla_prod, r, N)$estimate
  x[4,ct] = nu
  x[5:6,ct] = nlm(pbla_gsem, c(2,1), pbla_prod, ru, N)$estimate
  rp = rpseudor(ru, nu / n)
  x[7,ct] = length(rp)
  x[8:9,ct] = nlm(pbla_gsem, c(2,1), pbla_prod, rp, N)$estimate
  Nu = ceiling(N * nu / n)
  x[10:11,ct] = nlm(pbla_gsem, c(2,1), pbla_prod, ru, Nu)$estimate
  
  if(!(ct %% U)){
    print(ct)
  }
}

y = x[,order(x[2,]/x[3,])] # sorted by R0
z = x[,order(x[4,]/x[1,])] # sorted by reporting

# plotting

acol = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
acol = c(acol[6], acol[7], acol[4], acol[8], acol[1])

# bias in R0

par(mfrow=c(1,1))
plot(y[2,]/y[3,], y[2,]/y[3,] - y[5,]/y[6,], 
     pch = 19,cex = 0.5, col = acol[2],
     ylim = c(-0.1,1.5),
     xlab = expression("R"[0]), ylab = "Bias")
lines(y[2,]/y[3,], y[2,]/y[3,] - y[8,]/y[9,], col = acol[1], lwd = 2, lty = 1)
lines(y[2,]/y[3,], y[2,]/y[3,] - y[10,]/y[11,], col = acol[3], lwd = 1, lty = 1)
legend("topleft",
       legend = c("Underreported","Pseudo", "Scaled"),
       pch = rep(19,3),
       lty = rep(1,3),
       col = acol[c(2,1,3)])

# ratio difference in gamma and beta

# gamma

par(mfrow = c(2,1))
plot(z[4,]/z[1,], z[2,] / z[8,], 
     col = acol[1], type = 'l', lty = 3, 
     ylab = expression(beta), xlab = "Percent Reported",
     ylim = c(0.5,6))
lines(z[4,]/z[1,], z[2,]/z[5,], col = acol[2], lty = 3)
lines(z[4,]/z[1,], z[2,]/z[10,], col = acol[3], lty = 3)
abline(h=1, lty = 2)

d = z[2,]/z[5,]
p = z[4,]/z[1,]
bfl = loess(d ~ p)
bfl = predict(bfl, p)
lines(p, bfl, type = 'l', col = acol[2], lwd = 3)

d = z[2,]/z[8,]
p = z[4,]/z[1,]
bfl = loess(d ~ p)
bfl = predict(bfl, p)
lines(p, bfl, col = acol[1], lwd = 3)

d = z[2,]/z[10,]
p = z[4,]/z[1,]
bfl = loess(d ~ p)
bfl = predict(bfl, p)
lines(p, bfl, col = acol[3], lwd = 3)

legend("topright",
       c("Underreported", "Pseudo", "Scaled"),
       lwd = c(3,3,3),
       lty = c(1,1,1),
       col = acol[c(2,1,3)])

# gamma

plot(z[4,]/z[1,], z[3,]/z[9,], 
     col = acol[1], type = 'l', lty = 3,
     ylab = expression(gamma), xlab = "Percent Reported",
     ylim = c(0.5,3))
lines(z[4,]/z[1,], z[3,]/z[6,], col = acol[2], lty = 3)
lines(z[4,]/z[1,], z[3,]/z[11,], col = acol[3], lty = 3)
abline(h=1, lty = 2)

d = z[3,]/z[6,]
p = z[4,]/z[1,]
bfl = loess(d ~ p)
bfl = predict(bfl, p)
lines(p, bfl, type = 'l', col = acol[2], lwd = 3)

d = z[3,]/z[9,]
p = z[4,]/z[1,]
bfl = loess(d ~ p)
bfl = predict(bfl, p)
lines(p, bfl, col = acol[1], lwd = 3)

d = z[3,]/z[11,]
p = z[4,]/z[1,]
bfl = loess(d ~ p)
bfl = predict(bfl, p)
lines(p, bfl, col = acol[3], lwd = 3)
