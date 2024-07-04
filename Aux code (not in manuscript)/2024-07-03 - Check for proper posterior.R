# PRELIMINARIES ----------------------------------------------------

#rm(list=ls())

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
# renv::restore()

toLoad = c("crayon",
           "dplyr",
           "foreach",
           "doParallel",
           "metafor",
           "robumeta",
           "data.table",
           "purrr",
           "metRology",
           "fansi",
           "MetaUtility",
           "ICC",
           "cfdecomp",
           "tidyr",
           "tibble",
           "testthat",
           "rstan", # note: to reinstall this one, need to use high-mem session
           "optimx",
           "weightr",
           "phacking",
           "here",
           "stringr",
           "bayesmeta")  

# to install everything
# lapply(toLoad, install.packages)

lapply( toLoad,
        require,
        character.only = TRUE)



# prevent masking
select = dplyr::select

# run this only if you want to update the R environment specs
# setwd(here())
# renv::snapshot()

# ~~ User-specified global vars -------------------------
# no sci notation
options(scipen=999)


# ~~ Set directories -------------------------

# need helper fns from simulation study
code.dir = here()

setwd(code.dir)
source("helper_JTE.R")  # for lprior(), etc.


# DATASET ----------------------------------------------------

# this is one subset from the applied example
d = data.table::data.table(yi = c(-1.13943428318836, -1.10866262452161, -0.22314355131421),
                           vi = c(1.33090760913438, 2.65828884428493, 0.0342839539605235),
                           group = rep("All-cause death (k = 3)", 3L),
                           sei = c(1.1536496908223, 1.63042597019458, 0.185159266472201) )



# use nlpost function from helper_JTE!
# nlpost <- function(mu, tau, yi, sei) {
#   joint_nll <- get_nll(mu, tau, yi, sei) # negative log-likelihood
#   joint_lprior <- get_lprior(mu, tau, sei) # log-prior
#   joint_nll - joint_lprior # log-posterior
# }

# bm: make a wrapper fn in only mu, tau and do double integral over those two3 
ds = d[1,]
#ds=d
nlpost(mu = 0, tau = 0.01, yi = ds$yi, sei = ds$sei)

#** note: if tau = 0, then determinant is zero, so prior is infinity and so is posterior


ds = d[1,]


library(rmutil)

# the actual integral
int2(f = function(x, y) nlpost2(mu = x, tau = y, yi = ds$yi, sei = ds$sei),
          
          # lower bounds of integration (mu, tau):
          a = c(-Inf, -Inf),
          
          # upper bounds of integration:
          b = c(Inf, Inf))

# test with finite bounds
int2(f = function(x, y) nlpost2(mu = x, tau = y, yi = ds$yi, sei = ds$sei),
     
     # lower bounds of integration (mu, tau):
     a = c(0, 0),
     
     # upper bounds of integration:
     b = c(1, 1))

