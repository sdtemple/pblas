# PBLA Figure
# Seth Temple, sdtemple@uw.edu
# May 2021

# complementary to `pbla-exp.R`, `pbla-erl.R`, or `pbla-bakeoff.R`
# check array sizes
x = readRDS("some-epidemic.rds")
y = readRDS("another-epidemic.rds")

# colors
acol = c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
acol = c(acol[6], acol[7], acol[4], acol[8], acol[1])

par(mfrow=c(2,2))
xa = 0 # bounds
ya = 0

# Small -------------------------------------------------------------------

# beta

b = 0.3
yb = max(density(x[2,], from = 0)$y,
         density(x[4,], from = 0)$y,
         density(x[6,], from = 0)$y) + 0.1
xb = b * 2

plot.new()
plot.window(c(xa, xb), c(ya, yb))
axis(1)
axis(2)
title("A")
title(xlab = expression(beta))
title(ylab = "Density")
box()

lines(density(x[2,], from = 0), col = acol[1], lwd = 2)
lines(density(x[4,], from = 0), col = acol[2], lwd = 2, lty = 2)
lines(density(x[6,], from = 0), col = acol[3], lwd = 2, lty = 3)
abline(v = b, col = 'black', lwd = 2)
legend("bottomright",
       c("Std","Prod","E+D"),
       col = acol[1:3],
       lwd = rep(2,2),
       lty = 1:3)

# gamma

g = 0.2
yb = max(density(x[3,], from = 0)$y,
         density(x[5,], from = 0)$y,
         density(x[7,], from = 0)$y) + 0.1
xb = g * 2

plot.new()
plot.window(c(xa, xb), c(ya, yb))
axis(1)
axis(2)
title(main='B')
title(xlab = expression(gamma))
title(ylab = "Density")
box()

lines(density(x[3,], from = 0), col = acol[1], lwd = 2)
lines(density(x[5,], from = 0), col = acol[2], lwd = 2, lty = 2)
lines(density(x[7,], from = 0), col = acol[3], lwd = 2, lty = 3)
abline(v = g, col = 'black', lwd = 2)

# Large --------------------------------------------------------------------

# beta

b = 3
yb = max(density(y[2,], from = 0)$y,
         density(y[4,], from = 0)$y,
         density(y[6,], from = 0)$y) + 0.1
xb = b * 2

plot.new()
plot.window(c(xa, xb), c(ya, yb))
axis(1)
axis(2)
title("C")
title(xlab = expression(beta))
title(ylab = "Density")
box()

lines(density(y[2,], from = 0), col = acol[1], lwd = 2)
lines(density(y[4,], from = 0), col = acol[2], lwd = 2, lty = 2)
lines(density(y[6,], from = 0), col = acol[3], lwd = 2, lty = 3)
abline(v = b, col = 'black', lwd = 2)

# gamma

g = 2
yb = max(density(y[3,], from = 0)$y,
         density(y[5,], from = 0)$y,
         density(y[7,], from = 0)$y) + 0.1
xb = g * 2

plot.new()
plot.window(c(xa, xb), c(ya, yb))
axis(1)
axis(2)
title(main='D')
title(xlab = expression(gamma))
title(ylab = "Density")
box()

lines(density(y[3,], from = 0), col = acol[1], lwd = 2)
lines(density(y[5,], from = 0), col = acol[2], lwd = 2, lty = 2)
lines(density(y[7,], from = 0), col = acol[3], lwd = 2, lty = 3)
abline(v = g, col = 'black', lwd = 2)
