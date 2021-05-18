# PBLA Simulation Study (n/N)
# Seth Temple, sdtemple@uw.edu

library(pblas)

# Set Parameter Values ----------------------------------------------------------

set.seed(5132021)

# epidemic
N = 200
A = 5 # at least this many infecteds
beta = 1.5
gamma = 1

# Simulation Study ----------------------------------------------------------

# initialize
ct = 0
K = 2000 # simulation size
U = 50 # print feedback
storage = array(NA, dim = c(8,K))

# run
while(ct < K){
  epi = rgsem(beta, gamma, N)
  r = (epi[,2])[is.finite(epi[,2])]
  r = sort(r)
  while((length(r) < A)){
    epi = rgsem(beta, gamma, N)
    r = (epi[,2])[is.finite(epi[,2])]
    r = sort(r)
  }
  ct = ct + 1
  
  storage[1,ct] = length(r)
  storage[2,ct] = N
  storage[3:4,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_std_gsem, r=r, N=N)$estimate
  storage[5:6,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_prod, r=r, N=N)$estimate
  storage[7:8,ct] = nlm(pbla_gsem, c(1,1), pbla=pbla_weak, r=r, N=N)$estimate
  
  
  if(!(ct %% U)){
    print(ct)
  }
}

# saving
saveRDS(storage, "pbla-nN.rds")

# Inference by Decile -----------------------------------------------------

storage = readRDS("pbla-nN.rds")

N = storage[2,1]
deciles = N * seq(.1, 1, by=.1)
deciles[10] = deciles[10] + 1

xa = ya = 0
acol = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
acol = c(acol[6], acol[7], acol[4], acol[8], acol[1])

# all eight
par(mfrow=c(3,2))
for(i in 1:6){
  
  # subset by deciles
  x = (storage[1,] >= deciles[i]) & (storage[1,] < deciles[i+1])
  y = storage[,x]
  print(dim(y)[2])
  
  # betas
  yb = max(density(y[3,], from = 0)$y,
           density(y[5,], from = 0)$y) + 0.1
  xb = beta * 2
  
  plot.new()
  plot.window(c(xa, xb), c(ya, yb))
  axis(1)
  axis(2)
  title(main=paste('Percent Infected:', round(deciles[i] / N, 2), '-', round(deciles[i+1] / N, 2)))
  title(xlab = expression(beta))
  title(ylab = "Density")
  box()
  
  lines(density(y[5,], from = 0), col = acol[2], lwd = 2)
  lines(density(y[3,], from = 0), col = acol[1], lwd = 2) # std over
  abline(v = beta)
  if((i %in% c(1,5))){
    legend("topright",
           c("Std","Prod"),
           col = acol[1:2],
           lwd = rep(2,2))
  }
  
  # gammas
  yb = max(density(y[4,], from = 0)$y,
           density(y[6,], from = 0)$y) + 0.1
  xb = beta * 2
  
  plot.new()
  plot.window(c(xa, xb), c(ya, yb))
  axis(1)
  axis(2)
  title(xlab = expression(gamma))
  title(ylab = "Density")
  box()
  
  lines(density(y[6,], from = 0), col = acol[2], lwd = 2)
  lines(density(y[4,], from = 0), col = acol[1], lwd = 2) # std over
  abline(v = gamma)
}

par(mfrow=c(2,2))
dec = 3:6
d = dec[1]
xa = 0
xb = 3
for(i in dec){
  
  # subset by deciles
  x = (storage[1,] >= deciles[i]) & (storage[1,] < deciles[i+1])
  y = storage[,x]
  
  # betas
  yb = max(density(y[4,], from = 0)$y,
           density(y[6,], from = 0)$y) + 0.1
  
  plot.new()
  plot.window(c(xa, xb), c(ya, yb))
  axis(1)
  axis(2)
  title(main=paste('Percent Infected:', round(deciles[i] / N, 2), '-', round(deciles[i+1] / N, 2)))
  title(xlab = expression(beta))
  title(ylab = "Density")
  box()
  
  lines(density(y[5,], from = 0), col = acol[2], lwd = 2)
  lines(density(y[3,], from = 0), col = acol[1], lwd = 2) # std over
  #abline(v = beta)
  if(i == d){
    legend("topright",
           c("Std","Prod"),
           col = acol[1:2],
           lwd = rep(2,2))
  }
}

par(mfrow=c(2,2))
dec = 2:5
d = dec[1]
for(i in dec){
  
  # subset by deciles
  x = (storage[1,] >= deciles[i]) & (storage[1,] < deciles[i+1])
  y = storage[,x]
  
  # gammas
  yb = max(density(y[3,], from = 0)$y,
           density(y[5,], from = 0)$y,
           density(y[7,], from = 0)$y) + 0.1
  xb = beta * 2
  
  plot.new()
  plot.window(c(xa, xb), c(ya, yb))
  axis(1)
  axis(2)
  title(main=paste('Percent Infected:', round(deciles[i] / N, 2), '-', round(deciles[i+1] / N, 2)))
  title(xlab = expression(beta))
  title(ylab = "Density")
  box()
  
  lines(density(y[5,], from = 0), col = acol[2], lwd = 2)
  lines(density(y[3,], from = 0), col = acol[1], lwd = 2) # std over
  lines(density(y[7,], from = 0), col = acol[3], lwd = 2)
  #abline(v = gamma)
  if(i == d){
    legend("topright",
           c("Std","Prod","Weak"),
           col = acol[1:3],
           lwd = rep(2,3))
  }
}

par(mfrow=c(2,1))
cl = rcartocolor::carto_pal(7, "TealGrn")
dec = 1:7
d = dec[1]  
plot.new()
plot.window(c(xa, xb), c(0, 1))
axis(1)
axis(2)
title(xlab = expression(beta))
title(ylab = "Density")
box()
legend("topright",
       legend = paste(seq(d/10, d/10 + (length(dec)-1)/10, length.out = length(dec)),
                      '-', seq(d/10, d/10 + (length(dec)-1)/10, length.out = length(dec)) + 0.1),
       col = cl[dec],
       lwd = rep(2, length(dec)),
       title = "% Infected")
for(i in dec){
  
  # subset by deciles
  x = (storage[1,] >= deciles[i]) & (storage[1,] < deciles[i+1])
  y = storage[,x]
  
  # gammas
  lines(density(y[3,], from = 0), col = cl[i], lwd = 2) # std over
  abline(v = beta)
}
dec = 1:7
d = dec[1]  
plot.new()
plot.window(c(xa, xb), c(0, 1.5))
axis(1)
axis(2)
title(xlab = expression(gamma))
title(ylab = "Density")
box()
for(i in dec){
  
  # subset by deciles
  x = (storage[1,] >= deciles[i]) & (storage[1,] < deciles[i+1])
  y = storage[,x]
  
  # gammas
  lines(density(y[4,], from = 0), col = cl[i], lwd = 2) # std over
  abline(v = gamma)
}
