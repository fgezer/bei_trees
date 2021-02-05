library(spatstat)

load(url('https://github.com/fgezer/bei_trees/raw/main/bei_data.RData'))
source(url('https://github.com/fgezer/bei_trees/raw/main/lift_sim.R'))

L1 = lift2D(x = bei$x, y = bei$y, f = elev, nleft = 40, stage = NULL, 
            keepnbrs = F, Lmat = T, rw = c(0,1000,0,500), method = 'convex')

L2 = lift2D(x = bei$x, y = bei$y, f = elev, nleft = 40, stage = NULL, 
            keepnbrs = F, Lmat = T, rw = c(0,1000,0,500), method = 'rectangle')

L3 = lift2D(x = bei$x, y = bei$y, f = elev, nleft = 40, stage = NULL, 
            keepnbrs = F, Lmat = T, rw = c(0,1000,0,500), method = 'double')

L4 = lift2D(x = bei$x, y = bei$y, f = elev, nleft = 40, stage = NULL, 
            keepnbrs = F, Lmat = T, rw = c(0,1000,0,500), method = 'base')

L5 = lift2D(x = bei$x, y = bei$y, f = elev, nleft = 40, stage = NULL, 
            keepnbrs = F, Lmat = T, rw = c(0,1000,0,500), method = 'augm')


save(L1, L2, L3, L4, L5, file = 'out_bei.RData')
