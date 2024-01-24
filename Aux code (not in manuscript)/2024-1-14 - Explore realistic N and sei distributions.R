
# Goal: Try different N distributions and compare resulting sei distributions to those in SAPB-E
#  (see "2024-01-13 Empirical dist of SEs in SAPB-E")

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

# ~ Create my own N distributions  -------------------------------------------------

d = sim_meta(  k.pub = 100,
               t2a = 0.0001,
               Mu = 0,
               true.dist = "norm",
               p0 = 0.05,
               Ytype = "bin-OR",
               N.expr = "round( runif(n=1, min=40, max = 4000) )" )
hist(d$sei)


d = sim_meta(  k.pub = 100,
               t2a = 0.0001,
               Mu = 0,
               true.dist = "norm",
               p0 = 0.05,
               Ytype = "bin-OR",
               N.expr = "round( runif(n=1, min=40, max = 2000) )" )
hist(d$sei)


# CONTINUOUS Y  -------------------------------------------------

# try continuous Y
d = sim_meta(  k.pub = 100,
               t2a = 0.0001,
               Mu = 0.5,
               true.dist = "norm",
               Ytype = "cont-SMD",
               N.expr = "round( runif(n=1, min=2000, max = 4000) )"
               #N.expr = "round( runif(n=1, min=40, max = 4000) )"
)
hist(d$sei)


# try continuous Y
d = sim_meta(  k.pub = 100,
               t2a = 0.0001,
               Mu = 0.5,
               true.dist = "norm",
               Ytype = "cont-SMD",
               N.expr = "round( runif(n=1, min=400, max = 400) )"
               #N.expr = "round( runif(n=1, min=40, max = 4000) )"
)
hist(d$sei)









