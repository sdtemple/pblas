# PBLA Rabies in Central African Republic
# Seth Temple, sdtemple@uw.edu
# May 22, 2021

library(pblas)
library(outbreaks)
library(chron)
library(MASS)
library(ellipse)
library(stargazer)

# Data Prep ---------------------------------------------------------------

data("rabies_car_2003")
rab = rabies_car_2003$linelist
rab = rab[29:151,] # 2006-2012 epidemic

# time
rbegin = chron("2006-05-12", format=c(dates = "y-m-d")) 
rtimes = as.numeric(rab[,2] - rbegin)
rab$rt = rtimes

scalar = 28 # four weeks
r = rab$rt / scalar
plot(density(r), main = NA)

# Analysis ----------------------------------------------------------------

N = c(10000,25000,50000,100000)
k = 4 # which N
eta = c(.5, .2, .1)
table = matrix(0, 4, 6)
colnames(table) = c("eta","beta","gamma","R","Ra","Rb")

# 100% reporting
table[,1] = c(1, eta)
out = nlm(pbla_gsem, 
          c(2,1), 
          pbla=pbla_prod, 
          r=r, 
          N=N[k],
          hessian = T)
hessian = out$hessian
covariance = solve(hessian)
region = ellipse(covariance, centre = out$estimate, npoints = 1000)
table[1,2:3] = out$estimate
table[1,4] = table[1,2] / table[1,3]
x = seq(0.01, .5, length.out = 50)
y = seq(0.01, .5, length.out = 50)
z = matrix(0, 50, 50)
for(i in 1:50){
    print(i)
    for(j in 1:50){
        z[i,j] = - pbla_prod(r, x[i], y[j], N[k])
    }
}
filled.contour(x=x, y=y, z=z,
               xlab = expression(beta),
               ylab = expression(gamma),
               color.palette = terrain.colors,
               plot.axes = {
                   points(table[1,2], table[1,3], pch = 19)
                   axis(1, seq(0, .5, length.out = 6))
                   axis(2, seq(0, .5, length.out = 6))
                   points(region, pch = 19, cex = .5)
               })
title(main = "100% Reporting")
ratio = region[,1] / region[,2]
table[1,5] = min(ratio)
table[1,6] = max(ratio)

# 50% reporting
rp = rpseudor(r, eta[1])
rp = rp[rp >= min(r)]
rp = rp[rp <= max(r)]
#plot(density(rp), main = NA)
out = nlm(pbla_gsem, 
          c(2,1), 
          pbla=pbla_prod, 
          r=rp, 
          N=N[k],
          hessian = T)
hessian = out$hessian
covariance = solve(hessian)
region = ellipse(covariance, centre = out$estimate, npoints = 1000)
table[2,2:3] = out$estimate
table[2,4] = table[2,2] / table[2,3]
x = seq(.01, .5, length.out = 50)
y = seq(.01, .5, length.out = 50)
z = matrix(0, 50, 50)
for(i in 1:50){
    print(i)
    for(j in 1:50){
        z[i,j] = - pbla_prod(rp, x[i], y[j], N[k])
    }
}
filled.contour(x=x, y=y, z=z,
               xlab = expression(beta),
               ylab = expression(gamma),
               color.palette = terrain.colors,
               plot.axes = {
                   points(table[2,2], table[2,3], pch = 19)
                   axis(1, seq(0, .5, length.out = 6))
                   axis(2, seq(0, .5, length.out = 6))
                   points(region, pch = 19, cex = .5)
               })
title(main = "50% Reporting")
ratio = region[,1] / region[,2]
table[2,5] = min(ratio)
table[2,6] = max(ratio)

# 20% reporting
rp = rpseudor(r, eta[2])
rp = rp[rp >= min(r)]
rp = rp[rp <= max(r)]
#plot(density(rp), main = NA)
out = nlm(pbla_gsem, 
          c(2,1), 
          pbla=pbla_prod, 
          r=rp, 
          N=N[k],
          hessian=T)
hessian = out$hessian
covariance = solve(hessian)
region = ellipse(covariance, centre = out$estimate, npoints = 1000)
table[3,2:3] = out$estimate
table[3,4] = table[3,2] / table[3,3]
x = seq(.1, .6, length.out = 50)
y = seq(.1, .6, length.out = 50)
z = matrix(0, 50, 50)
for(i in 1:50){
    print(i)
    for(j in 1:50){
        z[i,j] = - pbla_prod(rp, x[i], y[j], N[k])
    }
}
filled.contour(x=x, y=y, z=z,
               xlab = expression(beta),
               ylab = expression(gamma),
               color.palette = terrain.colors,
               plot.axes = {
                   points(table[3,2], table[3,3], pch = 19)
                   axis(1, seq(.1, .6, length.out = 6))
                   axis(2, seq(.1, .6, length.out = 6))
                   points(region, pch = 19, cex = .5)
               })
title(main = "20% Reporting")
ratio = region[,1] / region[,2]
table[3,5] = min(ratio)
table[3,6] = max(ratio)

# 10% Reporting
rp = rpseudor(r, eta[3])
rp = rp[rp >= min(r)]
rp = rp[rp <= max(r)]
#plot(density(rp), main = NA)
out = nlm(pbla_gsem, 
          c(2,1), 
          pbla=pbla_prod, 
          r=rp, 
          N=N[k],
          hessian=T)
hessian = out$hessian
covariance = solve(hessian)
region = ellipse(covariance, centre = out$estimate, npoints = 1000)
table[4,2:3] = out$estimate
table[4,4] = table[4,2] / table[4,3]
x = seq(.1, .6, length.out = 50)
y = seq(.1, .6, length.out = 50)
z = matrix(0, 50, 50)
for(i in 1:50){
    print(i)
    for(j in 1:50){
        z[i,j] = - pbla_prod(rp, x[i], y[j], N[k])
    }
}
filled.contour(x=x, y=y, z=z,
               xlab = expression(beta),
               ylab = expression(gamma),
               color.palette = terrain.colors,
               plot.axes = {
                   points(table[4,2], table[4,3], pch = 19)
                   axis(1, seq(.2, .6, length.out = 5))
                   axis(2, seq(.2, .6, length.out = 5))
                   points(region, cex = .5, pch = 19)
               })
title(main = "10% Reporting")
ratio = region[,1] / region[,2]
table[4,5] = min(ratio)
table[4,6] = max(ratio)

# latex table
stargazer(table)
