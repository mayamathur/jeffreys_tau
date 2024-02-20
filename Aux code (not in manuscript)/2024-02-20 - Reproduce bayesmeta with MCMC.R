
# Goal: Make sure I know where the various measures of central tendency are calculated as I expect

# Conclusions: See project log

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
           "rma.exact")  

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
code.dir = here()

setwd(code.dir)
source("analyze_sims_helper_JTE.R")
source("helper_JTE.R")  # for lprior(), etc.
source("bayesmeta_edited.R")

# initialize stan model
source("init_stan_model_JTE.R")


# DATASET ----------------------------------------------------

# this is one subset from the applied example
d = data.table::data.table(yi = c(-1.13943428318836, -1.10866262452161, -0.22314355131421),
  vi = c(1.33090760913438, 2.65828884428493, 0.0342839539605235),
  group = rep("All-cause death (k = 3)", 3L),
  sei = c(1.1536496908223, 1.63042597019458, 0.185159266472201) )


# CHECK MYBAYESMETA VS. NEW PKG ----------------------------------------------------

# agrees exactly :D

# ~ mybayesmeta  -------------------------------------------------
m = mybayesmeta(y = d$yi,
                sigma = d$sei,
                tau.prior = "Jeffreys2",
                interval.type = "central")

# sanity check: plot posterior and prior
# prior is the dashed line
# plot(m, prior = TRUE)

# marginal (not joint) intervals
tau_ci = as.numeric( m$post.interval(tau.level=0.95) ) 
mu_ci = as.numeric( m$post.interval(mu.level=0.95) )

# this method doesn't do point estimation of inference for tau
list( stats = data.frame( 
  Mhat = m$MAP["joint", "mu"],
  Shat = m$MAP["joint", "tau"],
  MLo = mu_ci[1],
  MHi = mu_ci[2],
  SLo = tau_ci[1],
  SHi = tau_ci[2] ) )


# ~ new pkg  -------------------------------------------------

m = bayesmeta(y = d$yi,
                sigma = d$sei,
                tau.prior = "overallJeffreys",
                interval.type = "central")

# sanity check: plot posterior and prior
# prior is the dashed line
# plot(m, prior = TRUE)

# marginal (not joint) intervals
tau_ci = as.numeric( m$post.interval(tau.level=0.95) ) 
mu_ci = as.numeric( m$post.interval(mu.level=0.95) )

# this method doesn't do point estimation of inference for tau
list( stats = data.frame( 
  Mhat = m$MAP["joint", "mu"],
  Shat = m$MAP["joint", "tau"],
  MLo = mu_ci[1],
  MHi = mu_ci[2],
  SLo = tau_ci[1],
  SHi = tau_ci[2] ) )


# ~ MCMC  -------------------------------------------------

.yi = d$yi
.sei = d$sei
.stan.maxtreedepth = 25
.stan.adapt_delta = 0.995

init.fcn = function(o){ list(mu = 0,
                             tau = 0.1 ) }

options(mc.cores = parallel::detectCores())


post = sampling(stan.model,
                cores = 1,
                refresh = 0,
                data = list( k = length(.yi),
                             sei = .sei,
                             y = .yi ),
                
                control = list(max_treedepth = .stan.maxtreedepth,
                               adapt_delta = .stan.adapt_delta),
                
                init = init.fcn)


( postSumm = summary(post)$summary )
if (is.null(postSumm)) stop("In stan, postSumm is null")

# pull out best iterate to pass to MAP optimization later
ext = rstan::extract(post) # a vector of all post-WU iterates across all chains
best.ind = which.max(ext$log_post)  # single iterate with best log-posterior should be very close to MAP


# run optim for the pmode
Mhat.MaxLP = ext$mu[best.ind]
Shat.MaxLP = ext$tau[best.ind]
mle_fit <- mle_params(Mhat.MaxLP, Shat.MaxLP, d$yi, d$sei)
modes <- c(mle_fit@coef[["mu"]], mle_fit@coef[["tau"]])
optim_converged <- mle_fit@details$convergence == 0


# posterior means, posterior medians, modes, and max-LP iterate
# Mhat = c( postSumm["mu", "mean"],
#           median( rstan::extract(post, "mu")[[1]] ),
#           ext$mu[best.ind] )
# 
# Shat = c( postSumm["tau", "mean"],
#           median( rstan::extract(post, "tau")[[1]] ),
#           ext$tau[best.ind] )


### Reproduce each point estimate
# using m from bayesmeta above

# JOINT pmode - extremely similar
expect_equal(modes[1], m$MAP["joint", "mu"])
expect_equal(modes[2], m$MAP["joint", "tau"])

# how different are the marginal MAPs?
# pretty different for tau
m$MAP["joint", "mu"]; m$MAP["marginal", "mu"]
m$MAP["joint", "tau"]; m$MAP["marginal", "tau"]

# marginal pmean
# first check how stan calculates it in the first place
expect_equal( postSumm["mu", "mean"], mean( rstan::extract(post, "mu")[[1]] ) )
# now check against bayesmeta
expect_equal(postSumm["mu", "mean"], m$summary["mean", "mu"])
expect_equal(postSumm["tau", "mean"], m$summary["mean", "tau"])  #@QUITE DIFFERENT

# marginal pmed
expect_equal( median( rstan::extract(post, "mu")[[1]] ), m$summary["median", "mu"])
expect_equal( median( rstan::extract(post, "tau")[[1]] ), m$summary["median", "tau"])  #@AGAIN QUITE DIFFERENT



# ~ Visualize the joint and marginal posteriors  -------------------------------------------------

# https://mc-stan.org/bayesplot/articles/plotting-mcmc-draws.html
library(bayesplot)

mcmc_pairs(post, pars = c("mu", "tau"))

# doesn't look great:
mcmc_trace(post, pars = c("mu", "tau"), 
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_dens_overlay(post, pars = c("mu", "tau"))

# also doesn't look great:
mcmc_violin(post, pars = c("mu", "tau"), probs = c(0.1, 0.5, 0.9))


# https://cran.r-project.org/web/packages/bayesplot/vignettes/visual-mcmc-diagnostics.html
np_cp = nuts_params(post)
mcmc_parcoord(post, np = np_cp, pars = c("mu", "tau"))

