# PBLA Computation Time
# Seth Temple, sdtemple@uw.edu

library(pblas)

tbl = matrix(NA, nrow = 8, ncol = 6)
colnames(tbl) = c('n','N','Std','Prod','Weak','E+D')

b = 1.5
g = 1

# epidemic 1
x = readRDS('b1.5-g1-m1-n95-N200.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[1,1] = n
tbl[1,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[1,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[1,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[1,5] = t

t = system.time(pbla_ed_gsem(r, b, g, N))[3]
tbl[1,6] = t

# epidemic 2
x = readRDS('b1.5-g1-m1-n185-N500.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[2,1] = n
tbl[2,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[2,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[2,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[2,5] = t

t = system.time(pbla_ed_gsem(r, b, g, N))[3]
tbl[2,6] = t

# epidemic 3
x = readRDS('b1.5-g1-m1-n428-N1000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[3,1] = n
tbl[3,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[3,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[3,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[3,5] = t

t = system.time(pbla_ed_gsem(r, b, g, N))[3]
tbl[3,6] = t

# epidemic 4
x = readRDS('b1.5-g1-m1-n1483-N2500.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[4,1] = n
tbl[4,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[4,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[4,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[4,5] = t

t = system.time(pbla_ed_gsem(r, b, g, N))[3]
tbl[4,6] = t

# epidemic 5
x = readRDS('b1.5-g1-m1-n2830-N5000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[5,1] = n
tbl[5,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[5,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[5,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[5,5] = t

t = system.time(pbla_ed_gsem(r, b, g, N))[3]
tbl[5,6] = t

# epidemic 6
x = readRDS('b1.5-g1-m1-n5927-N10000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[6,1] = n
tbl[6,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[6,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[6,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[6,5] = t

# epidemic 7
x = readRDS('b1.5-g1-m1-n11819-N20000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[7,1] = n
tbl[7,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[7,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[7,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[7,5] = t

# epidemic 7
x = readRDS('b1.5-g1-m1-n29024-N50000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[8,1] = n
tbl[8,2] = N

t = system.time(pbla_std_gsem(r, b, g, N))[3]
tbl[8,3] = t

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[8,4] = t

t = system.time(pbla_weak(r, b, g, N))[3]
tbl[8,5] = t

# saving
saveRDS(tbl, 'pbla-time.rds')

# latex table
tbl = readRDS('pbla-time.rds')
library(stargazer)
stargazer(tbl)

# colors
acol = c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
acol = c(acol[6], acol[7], acol[4], acol[8], acol[1])

# plotting
par(mfrow=c(1,1))
plot(tbl[,1], tbl[,3],
     type = 'l',
     lwd = 2,
     ylab = "Time (s)",
     xlab = "n infecteds",
     col = acol[1])
lines(tbl[,1], tbl[,4],
      lwd = 2,
      lty = 2,
      col = acol[2])
lines(tbl[,1], tbl[,5],
      lwd = 2,
      lty = 3,
      col = acol[3])
lines(tbl[,1], tbl[,6],
      lwd = 2,
      lty = 4,
      col = acol[4])
legend("topleft",
       c("Std","Prod","Weak","E+D"),
       lwd = c(2,2,2,2),
       lty = c(1,2,3,4),
       col = acol[1:4])
