
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
           "stringr")  # note: to reinstall this one, need ml load jags

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

# control which results should be redone and/or overwritten
# but note that not all fns respect this setting
overwrite.res = TRUE

# should analyses be run from scratch?
# FALSE if just redoing plots
rerun.analyses = TRUE


# ~~ Set directories -------------------------

# need helper fns from simulation study
( code.dir = str_replace( string = here(),
                          pattern = "Applied example",
                          replacement = "Simulation study/Code \\(git\\)") )
# test it
setwd(code.dir)

data.dir = "/Users/mmathur/Dropbox/Personal computer/Independent studies/2024/JTE (Jeffreys estimation of tau)/Applied example"

# check that they're specified correctly
setwd(data.dir)

# below are the only absolute paths
# write results directly to directory containing TeX manuscript in Overleaf so stats can be piped directly into text
# this is an absolute path because it must live in Dropbox, outside the project directory, in order to sync with Overleaf
# to reproduce results, just set this to any directory on your local machine
# results will be written to a csv file in that location
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/JTE (Jeffreys tau estimation) Overleaf/R_objects/figures"
overleaf.dir.stats = "/Users/mmathur/Dropbox/Apps/Overleaf/JTE (Jeffreys tau estimation) Overleaf/R_objects"


setwd(code.dir)
source("analyze_sims_helper_JTE.R")
source("helper_JTE.R")  # for lprior(), etc.



# READ PREPPED DATA ----------------------------------------------------

setwd(data.dir)
d = fread("data_insomnia.csv")



.group = "TST; early follow-up (k = 4)"
d2 = d %>% filter( group == .group )


# MLE-profile: Look at optim methods  -------------------------------------------------

mle_params2 <- function(mu_start, tau_start, yi, sei) {
  nll_fun <- function(mu, tau) get_nll(mu, tau, yi, sei)
  stats4::mle(minuslogl = nll_fun,
              start = list(mu = mu_start, tau = tau_start),
              # using this method because in the single k=4 example below, 
              #  seems more likely than other methods to converge
              method = "L-BFGS-B")
}

my_mle = mle_params2(mu_start = 0,
                     tau_start = 0.1, 
                     yi = d2$yi,
                     sei = d2$sei)
my_mle
# warning that original solution hadn't converged:
cis = stats4::confint(my_mle) 


od = get_optimx_dataframe(.yi = d2$yi, 
                          .sei = d2$sei,
                          .mu.start = 0,
                          .tau.start = 0.1)


# which method won?
od2 = od %>% select( all_of( namesWith(dat = od, pattern = "nll") ) )
min(od2)



# MetaLik pkg: doesn't work -------------------------------------------------

library(metaLik)
mod2 = metaLik( yi ~ 1,
                data = d2,
                sigma2 = vi)

summ = summary(mod2)

confint( summary(mod2) )
inf = test.metaLik(mod2, param=1)

# try their example meta-analysis
data(education)
m <- metaLik(y~1, data=education, sigma2=sigma2)
summary(m)
confint(m)
test.metaLik(m, param=1)

profile(m)

# how did the Annals authors get CIs from this package??




# Metafor Q-profile -------------------------------------------------

# NOT the same as profile lkl

# can we just use metafor?
# no: it won't profile mu
# https://wviechtb.github.io/metafor/reference/profile.rma.html
# some underlying fns: https://github.com/cran/metafor/blob/master/R/profile.rma.uni.r
# uses Q-profile intervals
# from here (https://wviechtb.github.io/metafor/reference/confint.rma.html):
# For objects of class "rma.uni" obtained with the rma.uni function, a confidence interval for the amount of (residual) heterogeneity (i.e., 𝜏2
# ) can be obtained by setting random=TRUE (which is the default). The interval is obtained iteratively either via the Q-profile method or via the generalized Q-statistic method (Hartung and Knapp, 2005; Viechtbauer, 2007; Jackson, 2013; Jackson et al., 2014). The latter is automatically used when the model was fitted with method="GENQ" or method="GENQM", the former is used in all other cases.
#**this DOES change if original model is fit using REML
mod = rma( yi = d2$yi,
           vi = d2$vi,
           method = "ML",
           knha = FALSE )
( x = confint(mod) )
( tau.lb = x$random["tau","ci.lb"] )
( tau.ub = x$random["tau","ci.ub"] )
data.frame( tau.lb, tau.ub )


# **My own MLE-profile -------------------------------------------------


# my profile CI, as in doParallel
nll_fun <- function(mu, tau) get_nll(mu, tau, d2$yi, d2$sei)
my_mle = stats4::mle(minuslogl = nll_fun,
                     start = list(mu = 0, tau = 0.1),
                     method = "L-BFGS-B")

# this fn will complain if optimizer in my_mle hasn't found the true max
#  which is a good thing
cis = stats4::confint(my_mle)
cis
#??? why is the CI for tau symmetric around 0? doesn't make sense.

prof = stats4::profile(my_mle) 
plot( prof )


x = as.data.frame( attr(prof, "profile")$tau )
# okay, so the issue is just that it's trying to profile values of tau that are <0
plot( x$par.vals[, "tau"],
      x$z )



# ~ Exact -------------------------------------------------

library(rma.exact)
rma.exact(yi = d2$yi,
          vi = d2$vi)

# ~ bayesmeta  -------------------------------------------------


### Bayesmeta with Jeffreys prior on tau alone
library(bayesmeta)

# sanity check: should match MLE
#  because when mu.prior isn't specified, defaults to uniform
# yes, matches :)
m = bayesmeta(y = d2$yi,
              sigma = d2$sei,
              tau.prior = "uniform")

plot(m)


m = bayesmeta(y = d2$yi,
              sigma = d2$sei,
              tau.prior = "Jeffreys")
m
m$MAP["joint", "mu"]
as.numeric( m$post.interval(tau.level=0.95) )

m$post.interval(mu.level=0.95)
m$post.interval(mu.level=0.95, method = "central")


# ~ Editing bayesmeta  -------------------------------------------------

m1 = mybayesmeta(y = d2$yi,
              sigma = d2$sei,
              tau.prior = "Jeffreys",
              interval.type = "shortest")
m1

m2 = mybayesmeta(y = d2$yi,
                sigma = d2$sei,
                tau.prior = "Jeffreys2")
m2

# compare to my fn by using central interval
m2 = mybayesmeta(y = d2$yi,
                 sigma = d2$sei,
                 tau.prior = "Jeffreys2",
                 interval.type = "shortest")
m2
plot(m2)
# interesting: now the upper limit for tau is much higher

# comparing to my own (insomnia_forest_results), the "MAP marginal" from bayesmeta agrees very closely with my posterior modes :)
# however, the CIs are pretty different, presumably because bayesmeta gives marginal posterior summary?













