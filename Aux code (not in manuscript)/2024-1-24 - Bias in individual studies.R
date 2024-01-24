
# Goal: Confirm understanding of why some small-N, binary-Y scenarios have bias even in individual studies.

# PRELIMINARIES  -------------------------------------------------


# helper fns
code.dir = here()
setwd(code.dir)
source("helper_JTE.R")


# BINARY Y  -------------------------------------------------

# ~ Langan's N distributions  -------------------------------------------------

# this is reasonable, though less variation than seen in SAPB-E
d = sim_meta(  k.pub = 100,
               t2a = 0.0001,
               Mu = 0,
               true.dist = "norm",
               p0 = 0.05,
               Ytype = "bin-OR",
               N.expr = "round( runif(n=1, min=2000, max = 4000) )" )
hist(d$sei)

# fixed N = 400: reasonable shape, but low variation
d = sim_meta(  k.pub = 100,
               t2a = 0.0001,
               Mu = 0,
               true.dist = "norm",
               p0 = 0.05,
               Ytype = "bin-OR",
               N.expr = "round( runif(n=1, min=400, max = 400) )" )
hist(d$sei)

# fixed N = 40: NOT reasonable shape or variation
d = sim_meta(  k.pub = 100,
               t2a = 0.0001,
               Mu = 0,
               true.dist = "norm",
               p0 = 0.05,
               Ytype = "bin-OR",
               N.expr = "round( runif(n=1, min=40, max = 40) )" )
hist(d$sei)




