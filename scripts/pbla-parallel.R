# PBLA Computation Time
# Seth Temple, sdtemple@uw.edu

library(pblas)
library(parallel)
library(doParallel)
library(foreach)

cl = parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

tbl = matrix(NA, nrow = 8, ncol = 4)
colnames(tbl) = c('n','N','Serial','Parallel')

b = 1.5
g = 1

# epidemic 1
x = readRDS('b1.5-g1-m1-n95-N200.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[1,1] = n
tbl[1,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[1,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[1,4] = t

# epidemic 2
x = readRDS('b1.5-g1-m1-n185-N500.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[2,1] = n
tbl[2,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[2,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[2,4] = t

# epidemic 3
x = readRDS('b1.5-g1-m1-n428-N1000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[3,1] = n
tbl[3,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[3,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[3,4] = t

# epidemic 4
x = readRDS('b1.5-g1-m1-n1483-N2500.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[4,1] = n
tbl[4,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[4,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[4,4] = t

# epidemic 5
x = readRDS('b1.5-g1-m1-n2830-N5000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[5,1] = n
tbl[5,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[5,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[5,4] = t

# epidemic 6
x = readRDS('b1.5-g1-m1-n5927-N10000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[6,1] = n
tbl[6,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[6,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[6,4] = t


# epidemic 7
x = readRDS('b1.5-g1-m1-n11819-N20000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[7,1] = n
tbl[7,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[7,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[7,4] = t

# epidemic 7
x = readRDS('b1.5-g1-m1-n29024-N50000.rds')
N = nrow(x)
r = x[,2][is.finite(x[,2])]
n = length(r)

tbl[8,1] = n
tbl[8,2] = N

t = system.time(pbla_prod(r, b, g, N))[3]
tbl[8,3] = t

t = system.time(pbla_prod_parallel(r, b, g, N))[3]
tbl[8,4] = t

# saving
saveRDS(tbl, 'pbla-parallel.rds')

# latex table
tbl = readRDS('pbla-parallel.rds')
library(stargazer)
stargazer(tbl)

stopCluster(cl)
