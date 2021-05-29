# PBLA Ebola Virus (Contour)
# Seth Temple, sdtemple@uw.edu
# May 21, 2020

library(ellipse)

# Plotting ----------------------------------------------------------------

mle = readRDS("guinea-nlm.rds")
b0 = mle$estimate[1]
k0 = mle$estimate[2]

# temple, stockdale
bt = b0/1e3
kt = k0/1e3

# althaus (manually enter)
ba = .231
ka = .000712

some = readRDS("guinea-some.rds")
b = unique(some[,1])/1e3
k = unique(some[,2])/1e3
z = readRDS("guinea-contour.rds")

# zoom in (adjust accordingly)
#b = b[15:30]
#k = k[5:20]
#z = z[15:30,5:20]

# ellipse
hessian = mle$hessian
covariance = solve(hessian)
region = ellipse(covariance, centre = mle$estimate, npoints = 1000)
region = region / 1e3

# filled contour plot
filled.contour(b, k, z,
               xlab = expression(beta[0]),
               ylab = expression("k"[0]),
               color.palette = terrain.colors,
               plot.axes={
                 axis(1,seq(0,0.5,length.out=6)) # adjust accordingly
                 axis(2,seq(0,0.005, length.out=6)) # adjust accordingly
                 points(region, pch=19, cex=0.5)
                 points(bt, kt, pch=19, cex=2, col = "black")
                 points(ba, ka, pch=19, cex=2, col = "red")
                 })
title(main = "Guinea")

# study R0 ranges
min(region[,1]) * 5.61 # min
max(region[,1]) * 5.61 # max
