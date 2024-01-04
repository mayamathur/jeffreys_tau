
# NOTES ---------------------------------------------------------------

# "draw" always refers to within-study draw
# "study set" refers to all draws, observed and observed, for a given study


# ESTIMATION METHOD FNS ----------------------------------------------

# ~ Estimation Methods Structured for Use Inside run_method_safe --------

# Because these fns are run inside run_method_safe, the latter will handle editing rep.res
#  All of these fns should take get.CIs as an argument and return CIs as c(NA, NA) if not wanted
# Fns in this category need to return a dataframe with the below structure, although it's okay if they don't return all of these names since run_method_safe will handle that. Note that this is a LIST containing a dataframe called "stats", not just the dataframe; this allows easy extension in case you want to return other objects, like model objects.

# mu.start, sigma.start: start values for optimization
# as illustrated in a sanity check after nlpost_simple, this fn's MAPs agree with
#  using mle() directly on nlpost_Jeffreys
# confirmed that this agrees with weightr when usePrior = FALSE; see "2021-11-23 repurpose TNE code"
estimate_jeffreys_RTMA = function( yi, 
                                   sei,
                                   par2is = "Tt",
                                   Mu.start,
                                   par2.start,
                                   tcrit,  # vector of length k
                                   
                                   usePrior = TRUE,
                                   get.CIs,
                                   CI.method = "wald",
                                   
                                   run.optimx = FALSE ) {
  
  
  # #TEST ONLY
  # dpn = dp[ dp$affirm == FALSE, ]
  # yi = dpn$yi
  # sei = sqrt(dpn$vi)
  # par2is = "Tt"
  # Mu.start = p$Mu
  # par2.start = sqrt(p$t2a + p$t2w)
  # tcrit = qnorm(.975)
  # usePrior = FALSE
  # get.CIs = TRUE
  # CI.method = "wald"
  
  
  # ~~ Get MAP or MLE by Calling mle ----
  # IMPORTANT: This fn cannot be moved outside the scope of estimate_jeffreys
  #  because mle() is too dumb to allow extra args (e.g., x) to be passed,
  #  so it's forced to rely on global vars
  #  and that's a problem with a doParallel loop
  #  if this fn is outside estimate_jeffreys, different parallel iterations will use each other's global vars
  
  # ~~ Get MAP or MLE with Main Optimizer ----
  #  expects yi, sei, and tcrit to be global vars (wrt this inner fn)
  # second parameter could be Tt or T2t depending on par2is argument above
  nlpost_simple_RTMA = function(.Mu, .par2) {
    
    nlpost_jeffreys_RTMA( .pars = c(.Mu, .par2),
                          .par2is = par2is,
                          .yi = yi,
                          .sei = sei,
                          .tcrit = tcrit,
                          .usePrior = usePrior )
  }
  
  if ( usePrior == FALSE ) main.optimizer = "BFGS" else main.optimizer = "Nelder-Mead"
  
  mle.obj = mle( minuslogl = nlpost_simple_RTMA,
                 start = list( .Mu = Mu.start, .par2 = par2.start),
                 method = main.optimizer )
  
  
  # not actually MLEs, of course, but rather MAPs
  mles = as.numeric(coef(mle.obj))
  
  # THIS BEHAVES WELL
  if ( par2is == "Tt" ) {
    # need this structure for run_method_safe to understand
    MuHat = mles[1]
    TtHat = mles[2]
  }
  
  # from TNE
  # THIS BEHAVES BADLY
  if ( par2is == "T2t" ) {
    # need this structure for run_method_safe to understand
    MuHat = mles[1]
    TtHat = sqrt(mles[2])
  }
  
  # recode convergence more intuitively
  # optim uses "0" to mean successful convergence
  optim.converged = attr(mle.obj, "details")$convergence == 0 
  
  # ~~ Try Other Optimizers ----
  
  if ( run.optimx == TRUE ) {
    w = get_optimx_dataframe(.yi = yi,
                             .sei = sei,
                             .tcrit = tcrit,
                             .usePrior = usePrior,
                             .par2is = par2is,
                             .Mu.start = Mu.start,
                             .par2.start = par2.start )
  }
  
  
  # ~~ Inference ----
  # in case someone passes a set of params that aren't handled
  SEs = los = his = c(NA, NA)
  
  #profile.CI.error = NA
  
  if ( get.CIs == TRUE & CI.method == "wald" ) {
    # get Wald CI 
    # SEs for both parameters
    # these are from the observed Fisher info,
    #  as confirmed in "2021-8-9 Investigate inference.R"
    SEs = as.numeric( attr( summary(mle.obj), "coef" )[, "Std. Error" ] )
    
    # if needed, overwrite SE for second parameter so that it's 
    #  on Tt scale rather than T2t scale
    # delta method:
    # let g(y) = y^(1/2), where y=Tt
    #  so g'(y) = 0.5 * y^(-0.5)
    if ( par2is == "var" ) {
      TtHatSE = SEs[2] * 0.5*TtHat^(-0.5)
      SEs[2] = TtHatSE
    }
    
    los = mles - SEs * qnorm(0.975)
    his = mles + SEs * qnorm(0.975)
    
  } else if ( get.CIs == TRUE & CI.method != "wald" ) {
    warning("That CI.method isn't handled yet")
  }
  
  stats = data.frame( Mhat = MuHat, 
                      Shat = TtHat,
                      
                      MhatSE = SEs[1],
                      ShatSE = SEs[2],
                      
                      MLo = as.numeric(los[1]),
                      MHi = as.numeric(his[1]),
                      
                      SLo = as.numeric(los[2]),
                      SHi = as.numeric(his[2]),
                      
                      optim.converged = optim.converged )
  
  
  if ( run.optimx == TRUE ) stats = bind_cols(stats, w)
  
  return( list( stats = stats ) )
}



# .yi: published point estimates
# .sei: their SEs
# .tcrit: critical values on t or z scale for each study; can just use qnorm(.975) by default
# .Mu.start: optimizer starting value for meta-analysis Mu
# .Tt.start: optimizer starting value for meta-analysis tau
# .stan.adapt_delta: passed to rstan
# .stan.maxtreedepth: same
#  we should later use ellipsis to allow passing arbitrary args to rstan
estimate_jeffreys_mcmc_RTMA = function(.yi,
                                       .sei,
                                       .tcrit, 
                                       .Mu.start,
                                       .Tt.start,
                                       .stan.adapt_delta = 0.8,
                                       .stan.maxtreedepth = 10 ) {
  
  # stan.model (used later) is compiled OUTSIDE this fn in doParallel to avoid 
  #  issues with nodes competing with one another
  
  # handle scalar tcrit
  if ( length(.tcrit) < length(.yi) ) .tcrit = rep( .tcrit, length(.yi) )
  
  # prepare to capture warnings from Stan
  stan.warned = 0
  stan.warning = NA
  
  # set start values for sampler
  init.fcn = function(o){ list(mu = .Mu.start,
                               tau = .Tt.start ) }
  
  
  # like tryCatch, but captures warnings without stopping the function from
  #  returning its results
  withCallingHandlers({
    
    # necessary to prevent ReadRDS errors in which cores try to work with other cores' intermediate results
    # https://groups.google.com/g/stan-users/c/8snqQTTfWVs?pli=1
    options(mc.cores = parallel::detectCores())
    
    cat( paste("\n estimate_jeffreys_mcmc flag 2: about to call sampling") )
    
    post = sampling(stan.model,
                    cores = 1,
                    refresh = 0,
                    data = list( k = length(.yi),
                                 sei = .sei,
                                 tcrit = .tcrit,
                                 #2022-4-5: HANDLE AFFIRMATIVES FOR CSM CASE
                                 affirm = as.numeric( (.yi/.sei) > .tcrit ),
                                 y = .yi ),
                    
                    #iter = p$stan.iter,   
                    control = list(max_treedepth = p$stan.maxtreedepth,
                                   adapt_delta = p$stan.adapt_delta),
                    
                    init = init.fcn)
    
    
  }, warning = function(condition){
    stan.warned <<- 1
    stan.warning <<- condition$message
  } )
  
  cat( paste("\n estimate_jeffreys_mcmc flag 3: about to call postSumm") )
  postSumm = summary(post)$summary
  
  # pull out best iterate to pass to MAP optimization later
  ext = rstan::extract(post) # a vector of all post-WU iterates across all chains
  #best.ind = which.max(ext$lp__)  # single iterate with best log-posterior should be very close to MAP
  best.ind = which.max(ext$log_post)  # single iterate with best log-posterior should be very close to MAP
  
  
  # posterior means, posterior medians, and max-LP iterate
  Mhat = c( postSumm["mu", "mean"],
            median( rstan::extract(post, "mu")[[1]] ),
            ext$mu[best.ind] )
  
  Shat = c( postSumm["tau", "mean"],
            median( rstan::extract(post, "tau")[[1]] ),
            ext$tau[best.ind] )
  
  # sanity check
  expect_equal( Mhat[1], mean( rstan::extract(post, "mu")[[1]] ) )
  
  
  # SEs
  MhatSE = postSumm["mu", "se_mean"]
  ShatSE = postSumm["tau", "se_mean"]
  # how Stan estimates the SE: https://discourse.mc-stan.org/t/se-mean-in-print-stanfit/2869
  expect_equal( postSumm["mu", "sd"],
                sd( rstan::extract(post, "mu")[[1]] ) )
  expect_equal( MhatSE,
                postSumm["mu", "sd"] / sqrt( postSumm["mu", "n_eff"] ) )
  
  # CI limits
  S.CI = c( postSumm["tau", "2.5%"], postSumm["tau", "97.5%"] )
  M.CI = c( postSumm["mu", "2.5%"], postSumm["mu", "97.5%"] )
  # sanity check:
  myMhatCI = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], 0.025 ),
                            quantile( rstan::extract(post, "mu")[[1]], 0.975 ) ) )
  expect_equal(M.CI, myMhatCI)
  
  
  # the point estimates are length 2 (post means, then medians),
  #  but the inference is the same for each type of point estimate
  return( list( stats = data.frame( 
    
    Mhat = Mhat,
    Shat = Shat,
    
    MhatSE = MhatSE,
    ShatSE = ShatSE,
    
    # this will use same CI limits for all 3 pt estimates
    MLo = M.CI[1],
    MHi = M.CI[2],
    
    SLo = S.CI[1],
    SHi = S.CI[2],
    
    stan.warned = stan.warned,
    stan.warning = stan.warning,
    MhatRhat = postSumm["mu", "Rhat"],
    ShatRhat = postSumm["tau", "Rhat"] ),
    
    post = post,
    postSumm = postSumm ) )
  
}



# nicely report a metafor or robumeta object with optional suffix to denote which model
report_meta = function(.mod,
                       .mod.type = "rma",  # "rma" or "robu"
                       .suffix = "") {
  
  if ( !is.null(.mod) ) {
    
    
    if ( .mod.type == "rma" ) {
      tau.CI = tau_CI(.mod)
      .res = data.frame( .mod$b,
                         .mod$ci.lb,
                         .mod$ci.ub,
                         
                         sqrt(.mod$tau2),
                         tau.CI[1],
                         tau.CI[2] )
    } 
    
    
    if ( .mod.type == "robu" ) {
      
      .res = data.frame( .mod$b.r,
                         .mod$reg_table$CI.L,
                         .mod$reg_table$CI.U,
                         
                         sqrt(.mod$mod_info$tau.sq),
                         NA,
                         NA )
    } 
    
  } else {
    .res = data.frame( rep(NA, 6) )
  }
  
  
  names(.res) = paste( c("Mhat", "MLo", "MHi", "Shat", "SLo", "SHi"), .suffix, sep = "" )
  row.names(.res) = NULL
  
  return( list(stats = .res) )
}



# x: data vector
# doesn't handle the case CI.method = "profile" and par2is = "sd" (will just return NAs)
estimate_mle = function( x,
                         p,  # scenario parameters
                         par2is = "sd",
                         mu.start = 0,
                         sigma.start = 1,
                         get.CIs,  # this is kept separate from p for use inside boot fn (where we'll want to suppress getting CIs)
                         CI.method = "wald"
) {
  
  
  # ~~ Get MLE with Main Optimizer (NM) ----
  # fn needs to be formatted exactly like this (no additional args)
  #  in order for mle() to understand
  nll_simple = function(.mu, .sigma) {
    nll(.pars = c(.mu, .sigma),
        par2is = par2is,
        .x = x, .a = p$a, .b = p$b)
  }
  
  myMLE = mle( minuslogl = nll_simple,
               method = "BFGS",
               start = list( .mu=mu.start, .sigma=sigma.start) )
  
  mles = as.numeric( coef(myMLE) )
  
  
  # this parameterization behaves badly!
  if ( par2is == "sd" ) {
    # need this structure for run_method_safe to understand
    Mhat = mles[1]
    Vhat = mles[2]^2
    Shat = mles[2]
  }
  
  # this parameterization behaves well
  if ( par2is == "var" ) {
    # need this structure for run_method_safe to understand
    Mhat = mles[1]
    Vhat = mles[2]
    Shat = sqrt(mles[2])
  }
  
  mles = c(Mhat, Vhat, Shat)
  
  # recode convergence more intuitively
  # optim uses "0" to mean successful convergence
  optim.converged = attr(myMLE, "details")$convergence == 0
  
  
  # ~~ Try Other Optimizers ----
  w = get_optimx_dataframe(.method = "mle",
                           x = x,
                           p = p,  # scenario parameters
                           par2is = par2is,
                           mu.start = mu.start,
                           sigma.start = sigma.start)
  
  
  # ~~ Inference ----
  # in case someone passes a set of params that aren't handled
  SEs = los = his = c(NA, NA)
  
  profile.CI.error = NA
  
  if ( get.CIs == TRUE & CI.method == "wald" & par2is == "sd" ) {
    # get Wald CI 
    # SEs for both parameters
    # these are from the observed Fisher info,
    #  as confirmed in "2021-8-9 Investigate inference.R"
    SEs = as.numeric( attr( summary(myMLE), "coef" )[, "Std. Error" ] )
    
    # fill in VhatSE using delta method
    # let g(y) = y^2, where y=Shat
    VhatSE = SEs[2] * 2 * Shat
    
    # include Vhat
    SEs = c(SEs[1], VhatSE, SEs[2])
    
    los = mles - SEs * qnorm(0.975)
    his = mles + SEs * qnorm(0.975)
    
  } else if ( get.CIs == TRUE & CI.method == "wald" & par2is == "var" ) {
    # get Wald CI 
    # SEs for both parameters
    # these are from the observed Fisher info,
    #  as confirmed in "2021-8-9 Investigate inference.R"
    SEs = as.numeric( attr( summary(myMLE), "coef" )[, "Std. Error" ] )
    
    # fill in ShatSE using delta method
    # let g(y) = y^(1/2), where y=Shat
    #  so g'(y) = 0.5 * y^(-0.5)
    ShatSE = SEs[2] * 0.5*Vhat^(-0.5)
    
    # include Vhat
    SEs = c(SEs[1], SEs[2], ShatSE)
    
    los = mles - SEs * qnorm(0.975)
    his = mles + SEs * qnorm(0.975)
    
  } else if ( get.CIs == TRUE & CI.method == "profile" & par2is == "var" ) {
    
    tryCatch({
      CIs = confint(myMLE)
      
      los = c( as.numeric( CIs[1,1] ),
               as.numeric( CIs[2,1] ),
               sqrt( as.numeric( CIs[2,1] ) ) )
      his = c( as.numeric( CIs[1,2] ), 
               as.numeric( CIs[2,2] ),
               sqrt( as.numeric( CIs[2,2] ) ) )
      
      ( res.SEs = as.numeric( attr( summary(myMLE), "coef" )[, "Std. Error" ] ) )
      
      # fill in ShatSE using delta method
      # let g(y) = y^(1/2), where y=Shat
      #  so g'(y) = 0.5 * y^(-0.5)
      SEs = c(res.SEs[1], res.SEs[2], res.SEs[2] * 0.5*Vhat^(-0.5) )
      
    }, error = function(err) {
      SEs <<- los <<- his <<- rep(NA, 3)
      profile.CI.error <<- err$message
    })
    
  }   
  
  # ~~ Organize Results ----
  return( list( Mhat = Mhat, 
                Vhat = Vhat,
                Shat = Shat,
                
                MhatSE = SEs[1],
                VhatSE = SEs[2],
                ShatSE = SEs[3],
                
                M.CI = as.numeric( c(los[1], his[1]) ),
                V.CI = as.numeric( c(los[2], his[2]) ),
                S.CI = as.numeric( c(los[3], his[3]) ),
                optim.converged = optim.converged,
                profile.CI.error = profile.CI.error,
                
                optimx.dataframe = w
  ) )
}


# HELPERS FOR ABOVE ESTIMATION METHODS ----------------------------


# IMPORTANT NOTATION FOR THESE FNS:
# .Mu: mean of effect sizes, not Z-scores
# .T2t: total heterogeneity of effect sizes (=T2 + t2w in older notation, or t2a + t2w in newer notation)
# "T2t" is synonymous with "V" in parameters and results; "Tt" synonymous with "S" in params and results

# ~ Jeffreys prior fns ---------------

# **This fn can handle both nonaffirm and affirm results.
# agrees with weightr per "Repurpose TNE code.R"
# RTMA log-likelihood - now uses TNE version
# carefully structured for use with Deriv()
joint_nll_2 = function(.yi,
                       .sei,
                       .Mu,
                       .Tt = NULL,  # allow either parameterization
                       .T2t = NULL,
                       .tcrit = rep( qnorm(.975), length(.yi) ) ) {
  
  
  if ( is.null(.T2t) ) .T2t = .Tt^2
  if ( is.null(.Tt) ) .Tt = sqrt(.T2t)
  
  
  # as in TNE::nll()
  .dat = data.frame(yi = .yi,
                    sei = .sei,
                    crit = .tcrit,
                    affirm = (.yi/.sei) > .tcrit )
  
  .dat = .dat %>% rowwise() %>%
    mutate( term1 = dmvnorm(x = as.matrix(yi, nrow = 1),
                            mean = as.matrix(.Mu, nrow = 1),
                            # deal with different SEs by just incorporating them into the total variance
                            sigma = as.matrix(.T2t + sei^2, nrow=1),
                            log = TRUE),
            
            # this term applies if the study is nonaffirmative
            term2 = log( pmvnorm( lower = -99,
                                  upper = crit * sei,
                                  mean = .Mu,
                                  sigma = .T2t + sei^2 ) ),
            
            # this term applies if the study is affirmative
            # note that it's NOT 1-pmvnorm because it already handles the limits
            # e.g., pmvnorm(lower = 1.96, upper = Inf, mean = 0, sigma = 1) = 0.025
            #  (the upper tail)
            term3 = log( pmvnorm( lower = crit * sei,
                                  upper = 99,
                                  mean = .Mu,
                                  sigma = .T2t + sei^2 ) ),
            
            nll.i = ifelse( affirm == FALSE,
                            -term1 + term2,
                            -term1 + term3 ) )
  
  nll = sum(.dat$nll.i)
  
  # # sanity checks (from before allowing affirmatives)
  # term1.new = log( dnorm( x = .dat$yi[1],
  #                         mean = .Mu,
  #                         sd = sqrt(.T2t + .dat$sei[1]^2) ) )
  # 
  # term2.new = log( pnorm( q = .dat$crit[1] * .dat$sei[1],
  #                         mean = .Mu,
  #                         sd = sqrt(.T2t + .dat$sei[1]^2) ) ) 
  # 
  # expect_equal( as.numeric(.dat$term1[1]), as.numeric(term1.new), tol = 0.001 )
  # expect_equal( as.numeric(.dat$term2[1]), as.numeric(term2.new), tol = 0.001 )
  # 
  # # another sanity check
  # library(truncnorm)
  # nll.new = -log( dtruncnorm( x = .dat$yi[1],
  #                             mean = .Mu,
  #                             sd = sqrt(.T2t + .dat$sei[1]^2),
  #                             a = -99,
  #                             b = .dat$sei[1] * .dat$crit[1] ) )
  # 
  # expect_equal( nll.new, as.numeric(.dat$nll.i[1]) , tol = 0.001)
  
  # return it
  return(nll)
  
  # from TNE's nll:
  # term1 = dnorm(x = .x,
  #               mean = .mu,
  #               sd = .sigma,  
  #               log = TRUE)
  # 
  # term2 = length(.x) * log( pmvnorm(lower = .a,
  #                                   upper = .b,
  #                                   mean = .mu,
  #                                   # note use of sigma^2 here because of pmvnorm's different parameterization:
  #                                   sigma = .sigma^2 ) ) 
}

# ### Sanity check #2022-3-28b:
# # check against theory in paper
# d = sim_meta_2( Nmax = 1,
#                 Mu = 0.4,
#                 t2a = 0.4,
#                 m = 20,
#                 t2w = 0.01,
#                 true.sei.expr = "runif(n=1, 0.03, 1.5)",
#                 hack = "affirm",
#                 rho = 0,
#                 
#                 k.pub.nonaffirm = 10,
#                 prob.hacked = 0.8,
#                 return.only.published = FALSE)
# 
# Mu = 0.2
# Tt = 1
# ll1 = joint_nll_2(.yi = d$yi,
#                   .sei = sqrt(d$vi),
#                   .Mu = Mu,
#                   .Tt = Tt)
# # paper notation
# d$Si = sqrt(d$vi + Tt^2)  
# d$Vi = d$Si^2
# d$ctildei = (qnorm(0.975)*sqrt(d$vi) - Mu)/d$Si
# 
# 
# # just the normal itself without truncation
# # https://www.statlect.com/fundamentals-of-statistics/normal-distribution-maximum-likelihood
# ll1 = sum( dnorm(x = d$yi,
#                  mean = Mu, 
#                  sd = d$Si,
#                  log = TRUE) )
# ll2 = sum( -0.5*log(2*pi*d$Vi) - (1/(2*d$Vi))*(d$yi - Mu)^2 )
# expect_equal(ll1, ll2)
# 
# # now check the full RTMA
# dn = d %>% filter( yi/sqrt(vi) < qnorm(0.975) )
# ll1 = sum( log( dtruncnorm(x = dn$yi,
#                           a = -Inf,
#                           b = qnorm(0.975) * sqrt(dn$vi),
#                           mean = Mu,
#                           sd = dn$Si) ) )
# ll2 = sum( -0.5*log(2*pi*dn$Vi) - 1/(2*dn$Vi)*(dn$yi - Mu)^2 - log( pnorm(dn$ctildei) ) )
# expect_equal(ll1, ll2)
# ### end sanity check



# verbatim from TNE for sanity checks only
E_fisher_TNE = function(.mu, .sigma, .n, .a, .b) {
  
  # doesn't handle vectors because of max and min below
  if ( length(.sigma) > 1 ) stop("This fn doesn't handle multiple observations")
  
  # prevent infinite cutpoints
  # if either cutpoint is infinite, there are numerical issues because the alpha*Z terms
  #  below are 0*Inf
  Za = max( -99, (.a - .mu) / .sigma )
  Zb = min( 99, (.b - .mu) / .sigma )
  
  alpha.a = dnorm(Za) / ( pnorm(Zb) - pnorm(Za) )
  alpha.b = dnorm(Zb) / ( pnorm(Zb) - pnorm(Za) )
  
  k11 = -(.n/.sigma^2) + (.n/.sigma^2)*( (alpha.b - alpha.a)^2 + (alpha.b*Zb - alpha.a*Za) )
  
  k12 = -( 2*.n*(alpha.a - alpha.b) / .sigma^2 ) +
    (.n/.sigma^2)*( alpha.a - alpha.b + alpha.b*Zb^2 - alpha.a*Za^2 +
                      (alpha.a - alpha.b)*(alpha.a*Za - alpha.b*Zb) )
  
  k22 = (.n/.sigma^2) - (3*.n*(1 + alpha.a*Za - alpha.b*Zb) / .sigma^2) +
    (.n/.sigma^2)*( Zb*alpha.b*(Zb^2 - 2) - Za*alpha.a*(Za^2 - 2) +
                      (alpha.b*Zb - alpha.a*Za)^2 )
  
  return( matrix( c(-k11, -k12, -k12, -k22),
                  nrow = 2,
                  byrow = TRUE ) )
}


# Fisher info when taking derivatives wrt mu and tau
# This fn ONLY handles nonaffirm results, but could easily be adapted to handle
#  affirms by changing tcrit.
# important: note that in this fn, critical value is on t/z scale, NOT raw scale
#  vs. in E_fisher_TNE, .b is on raw scale
E_fisher_RTMA = function( .sei, .Mu, .Tt, .tcrit = qnorm(0.975) ) {
  
  Efish.list = lapply( X = as.list(.sei),
                       FUN = function(.s) {
                         
                         # for this observation
                         sei = .s
                         mu = .Mu
                         tau = .Tt
                         tcrit = .tcrit  # currently assumed to be a scalar
                         if ( length(tcrit) > 1 ) tcrit = tcrit[1] #OBVIOUSLY NEEDS TO BE GENERALIZED
                         
                         fishinfo = matrix( NA, nrow = 2, ncol = 2 )
                         
                         # from body of R's get_D11_num:
                         e2 = sei^2 + tau^2
                         e3 = sqrt(e2)
                         e5 = sei * tcrit - mu
                         e6 = e5/e3
                         e7 = dnorm(e6, 0, 1)
                         # Stan version:
                         # e7 = exp( normal_lpdf(e6 | 0, 1) )
                         e8 = pnorm(e6)
                         #e8 = exp( normal_lcdf(e6 | 0, 1 ) )
                         kmm = -(1/e2 - (e5/(e2 * e8) + e7 * e3/(e8 * e3)^2) * e7/e3)
                         
                         # from body of R's get_D12_num:
                         e2 = sei^2 + tau^2
                         e3 = sqrt(e2)
                         e5 = sei * tcrit - mu
                         # e6 is scaled critical value:
                         e6 = e5/e3
                         e7 = pnorm(e6)
                         # e7 = exp( normal_lcdf(e6 | 0, 1 ) )
                         e8 = e2^2
                         e9 = dnorm(e6, 0, 1)
                         #e9 = exp( normal_lpdf(e6 | 0, 1) )
                         
                         # my own expectation of .yi - .mu:
                         expectation1 = -sqrt(sei^2 + tau^2) * e9/e7
                         kms = -(tau * (((e7/e3 - e5 * e9/e2)/(e7 * e3)^2 - e5^2/(e8 *
                                                                                    e7 * e3)) * e9 + 2 * ((expectation1)/e8)))
                         
                         
                         # from body of R's get_D22_num:
                         e1 = tau^2
                         e3 = sei^2 + e1
                         e5 = sei * tcrit - mu
                         e6 = sqrt(e3)
                         # e7 is scaled crit value:
                         e7 = e5/e6
                         e8 = pnorm(e7)
                         # e8 = exp( normal_lcdf(e7 | 0, 1 ) )
                         e9 = dnorm(e7, 0, 1)
                         # e9 = exp( normal_lpdf(e7 | 0, 1 ) )
                         e10 = e5 * e9
                         e11 = e8 * e6
                         e13 = e10/e11
                         # *replace this one with its expectation:
                         # e15 = (.yi - .mu)^2/e3
                         # expectation of (.yi - .mu)^2:
                         expectation2 = (sei^2 + tau^2)*(1 - e7 * e9/e8)
                         e15 = (expectation2)/e3
                         
                         kss = (e13 + e15 - (e1 * (e5 * ((e8/e6 - e10/e3)/e11^2 -
                                                           e5^2/(e3^2 * e8 * e6)) * e9 + 2 * ((e13 + 2 * e15 -
                                                                                                 1)/e3)) + 1))/e3
                         
                         
                         fishinfo[1,1] = -kmm
                         fishinfo[1,2] = -kms
                         fishinfo[2,1] = -kms
                         fishinfo[2,2] = -kss
                         
                         return(fishinfo)
                         
                         
                       })
  
  # add all the matrices entrywise
  # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
  Efish.all = Reduce('+', Efish.list) 
  
  # cat("\nFirst observation Efish:")
  # print(Efish.list[[1]])
  
  return(Efish.all)
}

# ### Sanity check: Should agree with E_fisher_TNE if all SEs equal
# n = 10
# se = 2
# Mu = 0.1
# Tt = 2
# zcrit = 1.96
# ( Efish1 = E_fisher_RTMA( .sei = rep(se, n),
#                           .Mu = Mu,
#                           .Tt = Tt,
#                           .tcrit = rep(zcrit, n) ) )
# 
# ( Efish2 = E_fisher_TNE( .mu = Mu,
#                          .sigma = sqrt(Tt^2 + se^2),
#                          .n = n,
#                          .a = -99,
#                          .b = zcrit*se ) )
# expect_equal(Efish1, Efish2)
# 
# # c.f.: different SEs but with same mean across studies
# E_fisher_RTMA( .sei = runif(n = n, min = se - 1, max = se + 1),
#                .Mu = Mu,
#                .Tt = Tt,
#                .tcrit = rep(zcrit, n) )



lprior = function(.sei, .Mu, .Tt, .tcrit) {
  Efish = E_fisher_RTMA( .sei = .sei, .Mu = .Mu, .Tt = .Tt, .tcrit = .tcrit )
  log( sqrt( det(Efish) ) )
}


# with usePrior = FALSE, agrees with weightr per "Repurpose TNE code.R"
# .pars: (.Mu, .Tt) or (.Mu, .Tt2)
nlpost_jeffreys_RTMA = function( .pars,
                                 .par2is = "Tt",  # "Tt" or "T2t"
                                 .yi,
                                 .sei,
                                 .tcrit = qnorm(.975),
                                 
                                 # if .usePrior = FALSE, will just be the MLE
                                 .usePrior = TRUE) {
  
  
  if ( .pars[2] < 0 ) return(.Machine$integer.max)
  
  # variance parameterization
  if (.par2is == "T2t") {
    Mu = .pars[1]
    Tt = sqrt(.pars[2])
  }
  
  
  if (.par2is == "Tt") {
    Mu = .pars[1]
    Tt = .pars[2]
  }
  
  
  
  # negative log-likelihood
  # joint_nll_2 uses the TNE version
  nll.value = joint_nll_2( .yi = .yi,
                           .sei = .sei,
                           .Mu = Mu,
                           .Tt = Tt,
                           .tcrit = .tcrit )
  
  # log-prior
  # lprior uses the TNE expected Fisher and then just sums over observations
  if ( .usePrior == TRUE ) {
    prior.value = lprior(.sei = .sei,
                         .Mu = Mu,
                         .Tt = Tt,
                         .tcrit = .tcrit)
  } else {
    prior.value = 0
  }
  
  # negative log-posterior
  nlp.value = sum(nll.value) - prior.value
  
  if ( is.infinite(nlp.value) | is.na(nlp.value) ) return(.Machine$integer.max)
  
  return(nlp.value)
}


# ~ Other Helpers ---------------

# taken from TNE 2022-2-26
get_optimx_dataframe = function( .yi,
                                 .sei,
                                 .tcrit,
                                 .usePrior,
                                 .par2is,
                                 .Mu.start,
                                 .par2.start ) {
  
  
  ox.methods <- c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf',
                  'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
  
  l = optimx( par = c(.Mu.start, .par2.start),
              fn = function(..pars) as.numeric( nlpost_jeffreys_RTMA( .pars = ..pars,
                                                                      .par2is = .par2is,
                                                                      .yi = .yi,
                                                                      .sei = .sei,
                                                                      .tcrit = .tcrit,
                                                                      .usePrior = .usePrior ) ),
              method = ox.methods )
  
  l$opt.method = row.names(l)
  
  # transform second parameter so it's always Shat instead of Vhat
  if ( .par2is == "T2t" ) { l$p2 = sqrt(l$p2) }
  
  l2 = l %>% select(opt.method, p1, p2, convcode, value, kkt1, kkt2) 
  
  l2 = l2 %>% rename( Mhat = p1, Shat = p2, nll = value )
  
  w = pivot_wider(l2, 
                  names_from = "opt.method",
                  values_from = c("Mhat", "Shat", "convcode", "nll", "kkt1", "kkt2"),
                  names_glue = "optimx.{opt.method}.{.value}")
  
  
  if ( length( l$p1[ l$convcode == 0 ] ) > 0 ){
    
    # only keep the ones that had values for Mhat, Shat (not ones that didn't even give a value)
    l = l[ !is.na(l$p1) & !is.na(l$p2), ]
    
    
    #**optimizers that converged AND
    # had a small gradient (kkt1) AND
    # had a positive-definite Hessian (kkt2)
    lc = l[ l$convcode == 0 & l$kkt1 == TRUE & l$kkt2 == TRUE, ]
    
    # index of optimizer with the best nll
    lc.winner.ind = which.min(lc$value)
    
    
    # Mhat.winner is the Mhat of the optimizer with the best nll, OF converged ones
    # catch case in which no optimizers converged
    if ( length(lc.winner.ind > 0) ) {
      Mhat.winner = lc$p1[lc.winner.ind]
      Shat.winner = lc$p2[lc.winner.ind]
    } else {
      Mhat.winner = Shat.winner = NA
    }
    
    
    # **note that this is the criterion for agreement
    l$agree.Mhat = abs(l$p1 - Mhat.winner) < 0.01
    l$agree.Shat = abs(l$p2 - Shat.winner) < 0.01
    
    # sanity check: look at differences of non-agreers from Mhat.winner
    #l$p1[ l$agree.Mhat == FALSE ] - Mhat.winner
    
    
    # get lc again now that we have the agreement indicator
    lc = l[ l$convcode == 0 & l$kkt1 == TRUE & l$kkt2 == TRUE, ]
    w$optimx.Nconvergers = nrow(lc)
    w$optimx.convergers = paste( lc$opt.method, collapse = " ")
    
    w$optimx.Mhat.winner = Mhat.winner
    w$optimx.Pagree.Mhat.winner = sum(l$agree.Mhat)/nrow(l)
    # number and proportion of optimizers that converged that agreed with mode:
    w$optimx.Nagree.of.convergers.Mhat.winner = sum(lc$agree.Mhat)
    w$optimx.Pagree.of.convergers.Mhat.winner = sum(lc$agree.Mhat)/nrow(lc)
    w$optimx.Mhat.agreers = paste( l$opt.method[ l$agree.Mhat == TRUE ], collapse = " ")
    w$optimx.Mhat.convergers.agreers = paste( lc$opt.method[ lc$agree.Mhat == TRUE ], collapse = " ")
    
    w$optimx.Shat.winner = Shat.winner
    w$optimx.Pagree.Shat.winner = sum(l$agree.Shat)/nrow(l)
    w$optimx.Nagree.of.convergers.Shat.winner = sum(lc$agree.Shat)
    w$optimx.Pagree.of.convergers.Shat.winner = sum(lc$agree.Shat)/nrow(lc)
    w$optimx.Shat.agreers = paste( l$opt.method[ l$agree.Shat == TRUE ], collapse = " ")
    w$optimx.Shat.convergers.agreers = paste( lc$opt.method[ lc$agree.Shat == TRUE ], collapse = " ")
    
  } else {
    w$optimx.Nconvergers = NA
    w$optimx.convergers = NA
    w$optimx.Mhat.winner = NA
    w$optimx.Pagree.Mhat.winner = NA
    w$optimx.Nagree.of.convergers.Mhat.winner = NA
    w$optimx.Pagree.of.convergers.Mhat.winner = NA
    w$optimx.Mhat.agreers = NA
    w$optimx.Mhat.convergers.agreers = NA
    
    w$optimx.Shat.winner = NA
    w$optimx.Nagree.of.convergers.Shat.winner = NA
    w$optimx.Pagree.of.convergers.Shat.winner = NA
    w$optimx.Pagree.Shat.winner = NA
    w$optimx.Shat.agreers = NA
    w$optimx.Shat.convergers.agreers = NA
  }
  
  return(w)
} 


# ~ P-Curve Helpers --------------------

# Notes below are from McShane's R script:
# Code obtained on 02/10/2016 from:
# http://www.p-curve.com/Supplement/Rcode_paper2/APPENDIX%20-%20Loss%20Function%20and%20Estimation.R
# The code has been modified as per the following:
#	(i) The function originally named loss is now named pcurve.loss
#	(ii) The function originally named plotloss is now named pcurve.opt and plotting has been commented out.
#   (iii) The example has been commented out.

#LOSS FUNCTION
pcurve.loss=function(t_obs,df_obs,d_est) {
  #SYNTAX:
  #1. t_obs is a vector with observed t-values, 
  #2. df_obs vector with degrees of freedom associated with each t-value
  #3. d_est is the effect size on which fitted p-curve is based and the measure of loss computed
  
  
  #1.Convert all ts to the same sign (for justification see Supplement 5) 
  t_obs=abs(t_obs)
  
  #2 Compute p-values
  p_obs=2*(1-pt(t_obs,df=df_obs))
  
  #3 Keep significant t-values and corresponding df.
  t.sig=subset(t_obs,p_obs<.05)
  df.sig=subset(df_obs,p_obs<.05)
  
  
  #4.Compute non-centrality parameter implied by d_est and df_obs
  #df+2 is total N. 
  #Becuase the noncentrality parameter for the student distribution is ncp=sqrt(n/2)*d, 
  #we add 2 to d.f. to get N,  divide by 2 to get n, and by 2 again for ncp, so -->df+2/4
  ncp_est=sqrt((df.sig+2)/4)*d_est                          
  
  #5.Find critical t-value for p=.05 (two-sided)
  #this is used below to compute power, it is a vector as different tests have different dfs 
  #and hence different critical values
  tc=qt(.975,df.sig)                     
  
  #4.Find power for ncp given tc, again, this is a vector of implied power, for ncp_est,  for each test
  #MM: this can be 0 if d_est is very small
  power_est=1-pt(tc,df.sig,ncp_est)        
  
  #5.Compute pp-values
  #5.1 First get the overall probability of a t>tobs, given ncp
  p_larger=pt(t.sig,df=df.sig,ncp=ncp_est)
  
  #5.2 Now, condition on p<.05
  ppr=(p_larger-(1-power_est))/power_est  #this is the pp-value for right-skew
  
  # MM EDIT TO HANDLE EXTREME D_EST:
  # if anything goes wrong here, it sets the KSD p-value to 1
  #  such that this value of d_est will never be chosen as the min by pcurve.opt
  tryCatch({
    KSD=ks.test(ppr,punif)$statistic        #this is the D statistic outputted by the KS test against uniform
    return(KSD) 
  }, error = function(err){
    KSD = 1  # set KSD p-value to max possible value
    names(KSD) = "D"  # create same return structure as ks.test
    return(KSD)
  })
  
  # ORIGINAL VERSION (breaks if d_est is too extreme):
  # #6. Compute the gap between the distribution of observed pp-values and a uniform distribution 0,1 
  # KSD=ks.test(ppr,punif)$statistic        #this is the D statistic outputted by the KS test against uniform
  # return(KSD) 
}

#Function 2: Estimate d and plot loss function 
pcurve.opt=function(t_obs,df_obs,dmin,dmax){
  
  #SYNTAX:
  #t_obs  : vector with observed t-values 
  #df_obs : vector with degrees of freedom associated with each t-value
  #dmin   : smallest  effect size to consider 
  #dnax   : largest   effect size to consider 
  #e.g., dmin=-1, dmax=1 would look for the best fitting effect size in the d>=-1 and d<=1 range
  
  #Results will be stored in these vectors, create them first
  loss.all=c()
  di=c()
  
  #Compute loss for effect sizes between d=c(dmin,dmax) in steps of .01    
  for (i in 0:((dmax-dmin)*100))
  {
    d=dmin+i/100                   #effect size being considered
    di=c(di,d)                     #add it to the vector (kind of silly, but kept for symmetry)
    options(warn=-1)               #turn off warning becuase R does not like its own pt() function!
    loss.all=c(loss.all,pcurve.loss(df_obs=df_obs,t_obs=t_obs,d_est=d))
    #apply loss function so that effect size, store result
    options(warn=0)                #turn warnings back on
  }
  
  #find the effect leading to smallest loss in that set, that becomes the starting point in the optimize command
  imin=match(min(loss.all),loss.all)       #which i tested effect size lead to the overall minimum?
  dstart=dmin+imin/100                     #convert that i into a d.
  
  #optimize around the global minimum
  dhat=optimize(pcurve.loss,c(dstart-.1,dstart+.1), df_obs=df_obs,t_obs=t_obs)
  options(warn=-0)
  
  #Plot results
  #plot(di,loss.all,xlab="Effect size\nCohen-d", ylab="Loss (D stat in KS test)",ylim=c(0,1), main="How well does each effect size fit? (lower is better)")  
  ###points(dhat$minimum,dhat$objective,pch=19,col="red",cex=2)
  ###text(dhat$minimum,dhat$objective-.08,paste0("p-curve's estimate of effect size:\nd=",round(dhat$minimum,3)),col="red")
  return(dhat$minimum)
}


#Example
###t_obs= c(1.7,  2.8,  -3.1,  2.4)  #the first is ignored, the third converted to a positive number
###df_obs=c(44,   75,   125,  200)
###plotloss(t_obs=t_obs,df_obs=df_obs,dmin=-.5,dmax=2)

#MM examples
# pcurve.opt(t_obs=t_obs,df_obs=df_obs,dmin=-.5,dmax=2)
# pcurve.loss(t_obs=t_obs,df_obs=df_obs,d_est = 0)

# ANALYSIS FNS ---------------------------------------------------------------


# Notes from TNE:
# In order to catch errors from individual estimation methods safely and informatively,
#  in general the estimation method fns are structured st they can be run within the
#  fn run_method_safe, which automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error.

# ~~ Wrapper Fn to Safely Run a Method -------

# See note at the beginning of this script
#  this fn automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error

# Important: this fn works if method.fn() returns multiple rows
# BUT in that case, it assumes that the CIs are shared for all rows of that method

# expects global vars: all.errors, rep.res
# directly edits res via superassignment
run_method_safe = function( method.label,
                            method.fn,
                            .rep.res ) {
  
  cat( paste("\n run_method_safe flag 1: about to try running method", method.label) )
  
  
  tryCatch({
    
    method.output = method.fn()
    new.rows = method.output$stats
    
    if ( !exists("new.rows") ) {
      cat("\n\n**** Object new.rows didn't exist for method", method.label)
      cat("\nHere is method.output:\n")
      print(method.output)
    }
    
    cat( paste("\n run_method_safe flag 2: done calling method.fn() for", method.label) )
    
    error = NA
    
  }, error = function(err) {
    # needs to be superassignment because inside the "error" fn
    error <<- err$message
    
    # only need one variable in the blank dataframe since bind_rows below
    #  will fill in the rest
    new.rows <<- data.frame( method = method.label )
    
  })
  
  new.rows = new.rows %>% add_column( method = method.label, .before = 1 )
  new.rows$overall.error = error
  
  # optimx.dataframe is itself a df, so needs to be handled differently
  # if ( !is.null(optimx.dataframe) ) new.row = bind_cols(new.row, optimx.dataframe)
  
  if ( nrow(.rep.res) == 0 ) .rep.res = new.rows else .rep.res = bind_rows(.rep.res, new.rows)
  return(.rep.res) 
  
}

# example of how to call it when method.fn takes args
# all.errors = c()
# if  exists("rep.res") ) r("rep.re("rep.re
# run_method_safe( method = "mle",
#                  method.fn = function() estimate_mles(x = x, get.CIs = TRUE ) )

# #### Sanity checks
# # fake method for estimating the moments, but it breaks if x<0 or x>5
# crappy_method = function(x) {
#   if ( x > 0 & x < 5 ) return( list(Mhat = x+1,
#                                     Vhat = x-1,
#                                     M.CI = c(NA, NA),
#                                     V.CI = c(NA, NA) ) )
#   if ( x <= 0 ) stop("Fake error A generated by method.fn!")
#   if ( x >= 5 ) stop("Fake error B generated by method.fn!")
# }
# 
# all.errors = c()
# if( exists("rep.res") ) rm(rep.res)
# 
# run_method_safe( "mle",
#                  method.fn = function() { crappy_method(-1) } )
# 
# # no error on this one
# run_method_safe( "mle",
#                  method.fn = function() { crappy_method(4) } )
# 
# 
# # this one will have a different error
# # no error on this one
# run_method_safe( "mle", method.fn = function() { crappy_method(40) } )
# 
# expect_equal( all.errors, c( "mle: Fake error A generated by method.fn!",
#                             "mle: Fake error B generated by method.fn!" ) )
# 
# expect_equal( rep.res,
#               data.frame( method = rep("mle", 3),
#                           Mhat = c(NA, 5, NA),
#                           Vhat = c(NA, 3, NA),
#                           MLo = rep(NA, 3),
#                           MHi = rep(NA, 3),
#                           VLo = rep(NA, 3),
#                           VHi = rep(NA, 3) ) )
# #### end sanity checks


# ANALYSIS FNS: APPLIED EXAMPLES -----------------------------------------------


# ~ Fit Diagnostics -----------------------------------------------


# 2022-3-12
# fit diagnostics
# get CDF of (non-iid) marginal Z-scores (Zi.tilde)
#  given a fitted Shat
# .affirm: VECTOR with same length as x for affirm status
#  including the affirms is useful for 2PSM
yi_cdf = function(yi,
                  sei,
                  Mhat,
                  Shat) {
  
  
  affirm = (yi/sei) > qnorm(0.975)
  
  dat = data.frame( yi = yi,
                    sei = sei,
                    affirm = affirm,
                    cdfi = NA)
  
  if ( any(dat$affirm == FALSE) ) {
    dat$cdfi[ dat$affirm == FALSE ] = ptruncnorm(q = dat$yi[ dat$affirm == FALSE ],
                                                 a = -Inf,
                                                 b = qnorm(.975) * dat$sei[ dat$affirm == FALSE ],
                                                 mean = Mhat,
                                                 sd = sqrt(Shat^2 + sei[ dat$affirm == FALSE ]^2) )
  }
  
  if ( any(dat$affirm == TRUE) ) {
    dat$cdfi[ dat$affirm == TRUE ] = ptruncnorm(q = dat$yi[ dat$affirm == TRUE ],
                                                a = qnorm(.975) * dat$sei[ dat$affirm == TRUE ],
                                                b = Inf,
                                                mean = Mhat,
                                                sd = sqrt(Shat^2 + sei[ dat$affirm == TRUE ]^2))
  }
  
  return(dat)
  
}


# yi: published nonaffirmative estimates
# sei: their SEs
# Mhat: estimated Mu from RTMA
# Shat: estimated tau from RTMA
yi_qqplot = function(yi,
                     sei,
                     Mhat,
                     Shat){
  
  # get theoretical CDFs for each yi, given its affirm status
  cdf.dat = yi_cdf(yi = yi,
                   sei = sei,
                   Mhat = Mhat,
                   Shat = Shat)
  
  ecdf_fn = ecdf(yi)
  cdf.dat$ecdfi = ecdf_fn(yi)
  
  ggplot( data = cdf.dat,
          aes( x = cdfi,
               y = ecdfi) ) +
    geom_abline( slope = 1, 
                 intercept = 0,
                 color = "red") +
    geom_point( size = 2,
                alpha = 0.5 ) +
    xlab("Fitted CDF of point estimates") +
    ylab("Empirical CDF of point estimates") +
    theme_classic()
  
  
}

# ### sanity check on yi_cdf and yi_qqplot:
# Mhat=0.03033567
# Shat = 0.1010756
# setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/Lodder/Prepped data")
# dp = read.csv("lodder_prepped.csv")
# # published nonaffirmatives only
# dpn = dp[ dp$affirm == FALSE, ]
# kn= length(dpn$yi)
# # simulate data that conforms perfectly to fitted dist
# yin.fake = rtruncnorm( n = kn, 
#                        a = -Inf,
#                        b = qnorm(.975) * dpn$sei,
#                        mean = Mhat,
#                        sd = sqrt(Shat^2 + dpn$sei^2) )
# 
# 
# cdf.dat = yi_cdf(yi = yin.fake,
#                  sei = dpn$sei,
#                  Mhat = Mhat,
#                  Shat = Shat)
# 
# # sanity check
# cdfi.theory2 = ptruncnorm( yin.fake, 
#                            a = -Inf,
#                            b = qnorm(.975) * dpn$sei,
#                            mean = Mhat,
#                            sd = sqrt(Shat^2 + dpn$sei^2) )
# 
# expect_equal( cdf.dat$cdfi, cdfi.theory2 )
# 
# p = yi_qqplot(yi = yin.fake,
#               sei = dpn$sei,
#               Mhat = Mhat,
#               Shat = Shat)
# p
# ### end sanity check

# DATA SIMULATION ---------------------------------------------------------------


# runs a simple simulation to compare empirical moments to theoretical ones 
#  that I expected to match the empirical ones
# writes the dataset to .results.dir if it isn't NA (in which case needs to have a variable, .p$sim.name)
quick_sim = function(.p,
                     .results.dir = NA,
                     printRes = FALSE ) {
  
  # IMPORTANT: IF YOU ADD ARGS TO SIM_ONE_STUDY OR MAKE_ONE_DRAW, MUST ADD THEM 
  #  HERE OR ELSE QUICK_SIM WILL SILENTLY NOT PASS THEM
  # simulate a huge dataset, all finitely hacked
  d = sim_meta(Nmax = .p$Nmax,
               Mu = .p$Mu,
               T2 = .p$T2,
               m = .p$m,
               t2w = .p$t2w,
               se = .p$se,
               hack = .p$hack,
               return.only.published = FALSE,
               rho = .p$rho,
               
               k = .p$k,
               k.hacked = .p$k.hacked )
  
  
  # add in the parameters that aren't already in dataset
  shortParams = .p[ , !names(.p) %in% names(d) ]
  d = cbind( d, shortParams )
  
  # dataset of only published results
  dph = d[ d$Di == TRUE, ]
  
  Mu = unique(d$Mu)
  T2 = unique(d$T2)
  t2w = unique(d$t2w)
  m = unique(d$m)
  se = unique(d$se)
  
  library(msm)
  correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                             mean = c(Mu, se^2),
                             cov = matrix( c( T2 + t2w + se^2, 0, 0, var(d$vi) ),
                                           nrow = 2 ) )
  
  crit = unique(d$tcrit)
  
  # version in correct_meta_phack1 for nonaffirms:
  # for large m, var(d$vi) above is tiny, so sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ) is 
  # basically the same as correctedSE above
  # theoryExpTstat = extrunc(spec = "norm",
  #                          mean =.p$Mu /.p$se,
  #                          #doesn't use the delta-method thing
  #                          sd = sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ),
  #                          b = crit )
  
  
  res = data.frame( matrix( NA, nrow = 2, ncol = 2) )
  names(res) = c("affirms", "nonaffirms")
  row.names(res) = c("theoryExp", "empiricalExp")
  
  
  ## Expectation of affirmatives
  res[ "theoryExp", "affirms" ] = extrunc( spec = "norm",
                                           mean = Mu / se,
                                           sd = correctedSE,
                                           a = crit )
  
  res[ "empiricalExp", "affirms" ] = mean( dph$tstat[ dph$affirm == TRUE ] )
  
  
  ## Variance of affirmatives
  res[ "theoryVar", "affirms" ] = vartrunc( spec = "norm",
                                            mean = Mu / se,
                                            sd = correctedSE,
                                            a = crit )
  
  res[ "empiricalVar", "affirms" ] = var( dph$tstat[ dph$affirm == TRUE ] )
  
  
  # would be hard to look at affirms because of duplication within studies (messes up var)
  
  if ( printRes == TRUE ) print(res)
  
  library(Hmisc)
  returnList = llist(d, res, correctedSE)
  
  # save dataset
  if ( !is.na(.results.dir) & !is.na(.p$sim.name) ) {
    setwd(.results.dir)
    save( returnList, file=.p$sim.name )
  }
  
  
  return( returnList )
  
}



# changes from sim_meta:
# sei.expr
# k.pub.nonaffirm



# Simulate a meta-analysis in which some proportion of underlying studies (prior to publication)
#  are hacked, following various hacking mechanisms. Each study makes multiple draws (hypothesis tests)
#  until some stopping criterion based on Nmax and the hacking mechanism.

# - Nmax: max number of draws (hypothesis tests) that each hacked study can make before giving up
# - Mu: overall mean for meta-analysis
# - t2a: across-study heterogeneity (NOT total heterogeneity)
# Study parameters, assumed same for all studies:
#  - m: sample size
#  - t2w: within-study heterogeneity across draws
#  - true.sei.expr: quoted expression to evaluate to simulate a single study's standard error
#  - rho: autocorrelation of draws from a given study
#  - hack: mechanism of p-hacking (see sim_one_study_set for details)
# - k.pub.nonaffirm: number of published nonaffirmative studies desired in meta-analysis
#    (will simulate as many studies as needed to achieve that number)
# - prob.hacked: probability that an underlying study is hacked
sim_meta_2 = function(Nmax,  
                      Mu,  
                      t2a,  
                      
                      # study parameters, assumed same for all studies:
                      m,  # sample size for this study
                      t2w,  # within-study heterogeneity
                      true.sei.expr,  # TRUE SE string to evaluate
                      rho = 0,  # autocorrelation of muin's
                      
                      hack,  # mechanism of hacking for studies that DO hack (so not "no")
                      
                      k.pub.nonaffirm,  # number of published nonaffirmatives
                      prob.hacked,
                      
                      return.only.published = FALSE
) {
  
  
  # collect arguments
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k.pub.nonaffirm", "prob.hacked", "true.sei.expr")]
  
  
  if ( hack == "no" ) stop("hack should only be 'affirm' or 'signif' for this fn")
  
  k.pub.nonaffirm.achieved = 0
  i = 1
  
  while( k.pub.nonaffirm.achieved < k.pub.nonaffirm ) {
    
    is.hacked = rbinom(n = 1, size = 1, prob = prob.hacked)
    true.se = eval( parse( text = true.sei.expr) )
    
    # do we still need this??
    #if ( exists(".dat") ) startInd = max(.dat$study) + 1 else startInd = 1
    
    if ( is.hacked == 0 ) {
      # to generate unhacked studies, need to change argument "hack"
      .argsUnhacked = .args
      .argsUnhacked[ names(.args) == "hack" ] = "no"
      .argsUnhacked$se = true.se
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsUnhacked )
      
    } else if ( is.hacked == 1 ) {
      
      # for unhacked studies, no need to change argument "hack"
      .argsHacked = .args
      .argsHacked$se = true.se
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsHacked )
      
    }
    
    # add study ID
    newRows = newRows %>% add_column( .before = 1,
                                      study = i )
    
    # add study-draw ID
    study.draw = paste(newRows$study, 1:nrow(newRows), sep = "_")
    newRows = newRows %>% add_column( .after = 1,
                                      study.draw )
    
    if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    
    i = i + 1
    k.pub.nonaffirm.achieved = sum( .dat$affirm == FALSE & .dat$Di == 1 ) 
    
  }
  
  # add more info to dataset
  .dat$k.underlying = length(unique(.dat$study))
  .dat$k.nonaffirm.underlying = length( unique( .dat$study[ .dat$affirm == FALSE ] ) )
  
  if ( return.only.published == TRUE ) .dat = .dat[ .dat$Di == 1, ]
  
  return(.dat)
}

# d = sim_meta_2(  # test only
#   Nmax = 20,
#   Mu = 1,
#   t2a = 0.1,
#   m = 50,
#   t2w = .5,
#   true.sei.expr = "runif( n = 1, min = 0.5, max = 2 )",
#   hack = "affirm",
#   return.only.published = FALSE,
#   rho=0,
#   k.pub.nonaffirm = 30,
#   prob.hacked = 0.4 )
# 
# d$k.nonaffirm.underlying[1]
# d$k.underlying[1]
# table(d$Di, d$affirm)
# #KEEP/PULL IN THE SANITY CHECKS FROM THE VERSION BELOW


# KEEP THIS VERSION FOR REVERSE-COMPATIBILITY AND SANITY CHECKS
# *note that the number of reported, hacked studies might be less than k.hacked
#  if all Nmax draws are unsuccessful

# also note that for the unhacked study sets, the single published result could be 
# affirmative OR nonaffirmative

# simulate meta-analysis in which the hacking follows a mixture distribution:
# some studies are unhacked, in which case we always make Nmax draws and then report the last one (which is equivalent to only making 1 draw)
# and some studies are hacked, in which case we make UP TO Nmax draws and stop
# either when we get the first affirmative OR when we reach Nmax draws
sim_meta = function(Nmax,  # max draws to try
                    Mu,  # overall mean for meta-analysis
                    T2,  # across-study heterogeneity
                    
                    # study parameters, assumed same for all studies:
                    m,  # sample size for this study
                    t2w,  # within-study heterogeneity
                    se,  # TRUE SE
                    
                    rho = 0,  # autocorrelation of muin's
                    
                    hack,  # mechanism of hacking for studies that DO hack (so not "no")
                    
                    # args not passed to sim_one_study_set:
                    k,  # number of studies
                    k.hacked,  # number of hacked studies
                    
                    return.only.published = FALSE
) {
  
  # # test only
  # Nmax = 20
  # Mu = 0.1
  # T2 = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "affirm"
  # return.only.published = FALSE
  # rho=0
  # k = 30
  # k.hacked = 20
  
  # collect arguments
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k", "k.hacked")]
  
  
  if ( hack == "no" ) stop("hack should only be 'affirm' or 'signif' for this fn")
  
  k.unhacked = k - k.hacked
  
  
  ### Simulate the unhacked studies ###
  if ( k.unhacked > 0 ) {
    for ( i in 1:(k - k.hacked) ) {
      
      if ( i %% 50 == 0 ) cat("\nSimulating study #", i)
      
      # to generate unhacked studies, need to change argument "hack"
      .argsUnhacked = .args
      .argsUnhacked[ names(.args) == "hack" ] = "no"
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsUnhacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    }
  }
  
  ### Simulate hacked studies ###
  if ( k.hacked > 0 ) {
    if ( exists(".dat") ) startInd = max(.dat$study) + 1 else startInd = 1
    
    for ( i in startInd:(startInd + k.hacked - 1) ) {
      
      
      if ( i %% 50 == 0 ) cat("\nSimulating study #", i)
      
      # for unhacked studies, no need to change argument "hack"
      .argsHacked = .args
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsHacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    }
  }
  
  return(.dat)
  
}


# ### Example 1
# d = sim_meta(Nmax = 20,
#              Mu = 0.1,
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = 1,
#              hack = "affirm",
#              return.only.published = FALSE,
# 
#              k = 30,
#              k.hacked = 10
# 
# )
# 
# 
# 
# nrow(d)
# 
# # look at the published results only
# d %>% filter(Di == 1 ) %>%
#   group_by(hack) %>%
#   summarise( n(),
#              mean(affirm),
#              mean(mui),
#              mean(yi) )
# 
# ### Example 2: Affirm2 hacking
# d = sim_meta(Nmax = 5,
#              Mu = 0.1,
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = 1,
#              hack = "affirm2",
#              return.only.published = FALSE,
# 
#              k = 100,
#              k.hacked = 100
# 
# )
# 
# table(d$N)
# 
# # there should be some published nonaffirms
# d %>% filter(Di == 1) %>%
#   summarise( mean(affirm) )
# 
# # all are hacked, so every published nonaffirm should come from a study set
# #  that reached Nmax
# expect_equal( unique( d$N[ d$Di == 1 & d$affirm == FALSE] ),
#               5 )
# 
# 
# 
# nrow(d)
#
# ### Correlated draws
# # this is slow
# d = sim_meta(Nmax = 100,
#              Mu = -0.5,  # make it hard to get affirms
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = .5,
#              rho = 0.9,
#              hack = "affirm",
#              return.only.published = FALSE,
# 
#              k = 1000,
#              k.hacked = 500 )
# 
# # all results, even unpublished ones
# # primarily to check autocorrelation of muin's
# d %>% filter( !duplicated(study) ) %>%
#   group_by(hack) %>%
#   summarise( sum(!is.na(rhoEmp)),
#              sum(N > 1),
#              mean(rhoEmp, na.rm = TRUE) )


# ~ Simulate a single study ----------------- 

# Simulate study from potentially heterogeneous meta-analysis distribution;
#  within-study draws have their own heterogeneous within-study distribution

### Hacking types ###

# - "no": Makes exactly Nmax results and treats the last one as the reported one,
#     so the final result could be affirmative or nonaffirmative

# - "affirm" (worst-case hacking): Makes draws until the first affirmative is obtained
#    but if you reach Nmax, do NOT report any result at all. (Hack type "signif" is the 
#    same but favors all significant results.)

# - "affirm2": (NOT worst-case hacking): Similar to "affirm", but always reports the last draw,
#    even if it was nonaffirm (no file drawer)

# - "favor-best-affirm-wch" (worst-case hacking): Always makes Nmax draws. If you get any affirmatives,
#    publish the one with the lowest p-value. If you don't get any affirmatives, don't publish anything.

# NOTE: If you add args here, need to update quick_sim as well

# If Nmax is small, rhoEmp (empirical autocorrelation of muin's) will be smaller
#  than rho. That's okay because it reflects small-sample bias in autocorrelation
# estimate itself, not a problem with the simulation code
# For more about the small-sample bias: # https://www.jstor.org/stable/2332719?seq=1#metadata_info_tab_contents

sim_one_study_set = function(Nmax,  # max draws to try
                             Mu,  # overall mean for meta-analysis
                             t2a,  # across-study heterogeneity (NOT total heterogeneity)
                             m,  # sample size for this study
                             t2w,  # within-study heterogeneity
                             se,  # TRUE SE for this study
                             return.only.published = FALSE,
                             hack, # should this study set be hacked? ("no", "affirm","affirm2", "signif")
                             
                             # for correlated draws; see make_one_draw
                             rho = 0
) {  
  
  
  # # test only
  # Nmax = 20
  # Mu = 0.1
  # t2a = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "favor-best-affirm-wch"
  # rho=0
  
  # ~~ Mean for this study set ----
  # doesn't have t2w because that applies to results within this study set
  mui = Mu + rnorm(mean = 0,
                   sd = sqrt(t2a),
                   n = 1)
  
  # TRUE SD (not estimated)
  sd.y = se * sqrt(m)
  
  # collect all args from outer fn, including default ones
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  .args$mui = mui
  .args$sd.y = sd.y
  
  
  stop = FALSE  # indicator for whether to stop drawing results
  N = 0  # counts draws actually made
  
  # ~~ Draw until study reaches its stopping criterion ----
  # we use this loop whether there's hacking or not
  while( stop == FALSE & N < Nmax ) {
    
    if ( rho == 0 ) {
      # make uncorrelated draw
      newRow = do.call( make_one_draw, .args )
    } else {
      # make correlated draw
      if ( N == 0 ) .args$last.muin = NA  # on first draw, so there's no previous one
      if ( N > 0 ) .args$last.muin = d$muin[ nrow(d) ]
      newRow = do.call( make_one_draw, .args ) 
    }
    
    
    # number of draws made so far
    N = N + 1
    
    # add new draw to dataset
    if ( N == 1 ) d = newRow
    if ( N > 1 ) d = rbind( d, newRow )
    
    # check if it's time to stop drawing results
    if (hack == "signif") {
      stop = (newRow$pval < 0.05)
    } else if ( hack %in% c("affirm", "affirm2") ) {
      stop = (newRow$pval < 0.05 & newRow$yi > 0)
    } else if ( hack %in% c("no", "favor-best-affirm-wch") ) {
      # if this study set is unhacked, then stopping criterion
      #  is just whether we've reached Nmax draws
      # and for favor-best-affirm-wch, we always do Nmax draws so 
      #  we can pick the smallest p-value
      stop = (N == Nmax)
    } else {
      stop("No stopping criterion implemented for your chosen hack mechanism")
    }
    
  }  # end while-loop until N = Nmax or we succeed
  
  # record info in dataset
  d$N = N
  d$hack = hack
  
  # ~~ Empirical correlation of muin's ----
  #  but note this will be biased for rho in small samples (i.e., Nmax small)
  if ( nrow(d) > 1 ) {
    # get lag-1 autocorrelation
    d$rhoEmp = cor( d$muin[ 2:length(d$muin) ],
                    d$muin[ 1: ( length(d$muin) - 1 ) ] )
    
    # mostly for debugging; could remove later
    d$covEmp = cov( d$muin[ 2:length(d$muin) ],
                    d$muin[ 1: ( length(d$muin) - 1 ) ] )
    
  } else {
    d$rhoEmp = NA
    d$covEmp = NA
  }
  
  # convenience indicators for significance and affirmative status
  d$signif = d$pval < 0.05
  d$affirm = (d$pval < 0.05 & d$yi > 0)
  
  
  # ~~ Decide which draw to favor & publish ----
  # in the first 2 cases, Di=1 for only the last draw IF we got an affirm result
  #  but if we didn't, then it will always be 0
  if ( hack == "signif" ) d$Di = (d$signif == TRUE)
  if (hack == "affirm") d$Di = (d$affirm == TRUE)
  
  # if no hacking or affirmative hacking without file drawer,
  #   assume only LAST draw is published,
  #   which could be affirm or nonaffirm
  if ( hack %in% c("no", "affirm2") ) {
    d$Di = 0
    d$Di[ length(d$Di) ] = 1
  }
  
  # for favor-best-affirm-wch, favor the one with the 
  if ( hack %in% c("favor-best-affirm-wch") ) {
    d$Di = 0
    # if there was at least 1 affirm, publish it 
    if ( any(d$affirm == TRUE) ) {
      best.affirm.pval = min( d$pval[d$affirm == TRUE] )
      d$Di[ d$pval == best.affirm.pval & d$affirm == TRUE ] = 1
    }
    # ...otherwise don't publish any draw
    # sanity check:
    #View(d%>%select(Di,affirm,pval,yi))
  }
  
  if ( return.only.published == TRUE ) d = d[ d$Di == 1, ]
  
  return(d)
  
}


### example
# 
# d = sim_one_study_set(Nmax = 20,
#                       Mu = 1,
#                       t2a = 1,
#                       m = 50,
#                       t2w = .5,
#                       se = 1,
#                       hack = "favor-best-affirm-wch",
#                       return.only.published = FALSE)
# nrow(d)
# d



# ### example
# d = sim_one_study_set(Nmax = 5,
#                       Mu = 0.1,
#                       t2a = 0.1,
#                       m = 50,
#                       t2w = .5,
#                       se = 1,
#                       hack = "affirm2",
#                       return.only.published = FALSE)
# nrow(d)
# d


# ### sanity check by simulation
# for ( i in 1:2000 ) {
#   newRows = sim_one_study_set(Nmax = 20,
#                          Mu = 0.1,
#                          t2a = 0.1,
#                          m = 50,
#                          t2w = .1,
#                          se = 1,
#                          hack = "no",
#                          return.only.published = FALSE )
# 
#   if ( i == 1 ) .d = newRows else .d = rbind(.d, newRows)
# 
# }
# 
# # all studies
# # note that conditional on Di == 0, variance might be off because of repeated rows
# # but means should be correct
# .d %>% group_by(Di == 1) %>%
#   summarise(n(),
#             mean(mui),
#             mean(muin),
#             var(muin),
#             mean(yi) )
# # seems fine
# 
# ### check correlated draws: large Nmax
# # no hacking so that we can get a lot of draws
# d = sim_one_study_set(Nmax = 1000,
#                       Mu = 0,
#                       t2a = 0.1,
#                       m = 50,
#                       t2w = .5,
#                       se = 0.5,
#                       rho = 0.9,
#                       hack = "no",
#                       return.only.published = FALSE)
# nrow(d)
# # should match each other and should be close to rho
# table(d$rhoEmp)
# acf(d$muin, lag = 1)$acf
#
# ### check correlated draws: small Nmax
# #  wrong because of small-sample bias in autocorrelation (not my fn's fault)
# #  about the small-sample bias: # https://www.jstor.org/stable/2332719?seq=1#metadata_info_tab_contents
# # issue is that the sample variance and autocorrelation estimate are non-independent
# #  in small samples
# rhoEmp = c()
# covEmp = c()
# varEmp = c()
# for ( i in 1:250 ) {
#   d = sim_one_study_set(Nmax = 10,
#                         Mu = 0,
#                         t2a = 0.1,
#                         m = 50,
#                         t2w = .5,
#                         se = 0.5,
#                         rho = 0.9,
#                         hack = "no",
#                         return.only.published = FALSE)
#   
#   covEmp = c(covEmp, unique(d$covEmp))
#   rhoEmp = c(rhoEmp, unique(d$rhoEmp))
#   varEmp = c(varEmp, var(d$muin))
# }
# 
# # TOO LOW when Nmax is small! 
# mean(rhoEmp)  # when Nmax = 10: 0.46; when Nmax = 200: 0.88
# mean(covEmp)  # when Nmax = 10: 0.09; when Nmax = 200: 0.40
# mean(varEmp)  # when Nmax = 10: 0.15; when Nmax = 200: 0.47 (should equal t2w = 0.5)
# # sample variance is biased downward in small samples



# # ~ Sanity check  ---------------------------------------------------------------
# #  if Nmax -> Inf and we hack until affirmative, 
# #   published results should follow truncated t distribution
# # also need to set heterogeneity to 0?
# d = data.frame( matrix(nrow = 500, ncol = 1))
# Mu = 1
# t2a = 0.5
# t2w = 0.3
# se = 1
# 
# d = d %>% rowwise() %>%
#   mutate( sim_one_study_set( Nmax = Inf,
#                              Mu = Mu,
#                              t2a = t2a,
#                              m = 50,
#                              t2w = t2w,
#                              se = se,
#                              hack = "affirm",
#                              return.only.published = TRUE ) )
# 
# summary(d$N)
# 
# # calculate noncentrality parameter
# ncp = Mu / sqrt( t2a + t2w + se^2 )
# 
# # compare to truncated t
# qqtrunc(x = d$tstat,
#         spec = "t",
#         ncp = ncp,
#         df = 50-1,
#         # since I simulated iid studies here, the truncation cutoff is always the same
#         a = d$tcrit[1] )
# 
# # vs. actually drawing directly from truncated t
# x = rtrunc( n = 500,
#             spec = "t",
#             ncp = ncp,
#             df = 50-1,
#             a = d$tcrit[1])
# 
# plot( density(x) ) +
#   lines( density(d$tstat), col = "red")
# # looks great! :)



# ~ Draw one unbiased result within one study ------------------
# muin should be NA if either we want uncorrelated draws OR it's the first of a series of correlated draws
make_one_draw = function(mui,  # mean for this study set
                         t2w,
                         sd.y,  # TRUE SD
                         m,  # sample size
                         
                         # for making correlated draws
                         rho = 0,  # autocorrelation of muin's (not yi's)
                         last.muin = NA,  
                         ...) {
  
  
  # true mean for draw n (based on within-study heterogeneity)
  # either we want uncorrelated draws OR it's the first of a series of correlated draws
  if ( rho == 0 | is.na(last.muin) ) {
    muin = rnorm(mean = mui,
                 sd = sqrt(t2w),
                 n = 1)
  }
  
  # make correlated draw
  if ( rho != 0 & !is.na(last.muin) ) {
    # draw from BVN conditional, given muin from last draw
    # conditional moments given here:
    #  https://www.amherst.edu/system/files/media/1150/bivarnorm.PDF
    muin = rnorm(mean = mui + rho * (last.muin - mui),
                 sd = sqrt( t2w * (1 - rho^2) ),
                 n = 1)
  }
  
  
  # draw subject-level data from this study's population effect
  y = rnorm( mean = muin,
             sd = sd.y,
             n = m)
  
  # run a one-sample t-test
  test = t.test(y,
                alternative = "two.sided")
  
  pval = test$p.value
  tstat = test$statistic
  vi = test$stderr^2  # ESTIMATED variance
  
  
  # if (hack == "signif") success = (pval < 0.05)
  # if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
  
  return( data.frame(pval = pval,
                     tstat = tstat,
                     tcrit = qt(0.975, df = m-1),
                     mui = mui,
                     muin = muin,
                     yi = mean(y),
                     vi = vi,
                     viTrue = sd.y^2 / m,  # true variance; will equal p$se^2
                     m = m ) )
  #success = success,
  #N = Nmax ) )
}

# make_one_draw(mui = 0.1,
#               t2w = 0,
#               sd.y = 0.3,
#               m = 30 )
# 
# # rho = 1, so should get exactly the same muin again
# make_one_draw(mui = 0.1,
#               t2w = 0.1,
#               sd.y = 0.3,
#               m = 30,
#               rho = 1,
#               last.muin = 0.8)
# 
# # sanity check by simulation: uncorrelated draws
# for ( i in 1:5000 ) {
#   
#   # get last draw
#   last.muin = NA
#   if ( i > 1 ) last.muin = .d$muin[ nrow(.d) ]
#   
#   newRow = make_one_draw(mui = 0.1,
#                   t2w = 0.1,
#                   sd.y = 0.3,
#                   m = 30,
#                   rho = 0.9,
#                   last.muin = last.muin)
# 
#   if ( i == 1 ) .d = newRow else .d = rbind(.d, newRow)
# 
# }
# 
# .d %>% summarise( mean(mui),
#                   mean(muin),
#                   var(muin),
#                   mean(yi) )
# 
# 
# # sanity check by simulation: correlated draws
# for ( i in 1:5000 ) {
#   newRow = make_one_draw(mui = 0.1,
#                          t2w = 0.1,
#                          sd.y = 0.3,
#                          m = 30,
#                          rho = 1,
#                          last.muin = 0.8)
#   
#   if ( i == 1 ) .d = newRow else .d = rbind(.d, newRow)
#   
# }
# 
# .d %>% summarise(
#   mean(mui),
#   mean(muin),
#   var(muin),
#   mean(yi) )
# 
# # look at autocorrelation of muin's (should match rho above)
# acf(.d$muin, lag = 1)$acf
# 
# # look at autocorrelation of yi's (<rho because of SE>0)
# acf(.d$yi, lag = 1)$acf

# draw SE from Lodder distribution (with replacement)
#  this assumes that Lodder dataset has already been read in (by doParallel)
#  and there's global vector called "lodder.ses"
draw_lodder_se = function() {
  sample(x = lodder.ses,
         size = 1, 
         replace = TRUE)
}


# DATA SIMULATION: STEFAN ---------------------------------------------------------------


# reference:
# sim_one_study_set_stefan = function(strategy.stefan,
#                                     alternative.stefan,
#                                     
#                                     stringent.hack,
#                                     
#                                     return.only.published = FALSE,
#                                     is.hacked,  # separate from hack.type so that we can use same fn in each case for comparability
#                                     hack.type)


sim_meta_2_stefan = function(strategy.stefan,
                             alternative.stefan,
                             
                             stringent.hack,
                             prob.hacked,
                             hack.type,
                             k.pub.nonaffirm,
                             return.only.published = FALSE
) {
  

  # # test only
  # strategy.stefan = "firstsig"
  # alternative.stefan = "greater"
  # prob.hacked = 0.8
  # stringent.hack = TRUE
  # return.only.published = FALSE
  # hack.type = "DV"
  # k.pub.nonaffirm = 10  # global var in doParallel
  
  k.pub.nonaffirm.achieved = 0
  i = 1
  
  while( k.pub.nonaffirm.achieved < k.pub.nonaffirm ) {
    
    is.hacked = rbinom(n = 1, size = 1, prob = prob.hacked)
  
    
    newRow = sim_one_study_set_stefan(strategy.stefan = strategy.stefan,
                                      alternative.stefan = alternative.stefan,
                                      
                                      stringent.hack = stringent.hack,
                                      
                                      return.only.published = return.only.published,
                                      is.hacked = is.hacked,  
                                      hack.type = hack.type)
      
 
    
    # add study ID
    newRow = newRow %>% add_column( .before = 1,
                                      study = i )
  
    
    if ( i == 1 ) .dat = newRow else .dat = rbind( .dat, newRow )
    
    i = i + 1
    k.pub.nonaffirm.achieved = sum( .dat$affirm == FALSE & .dat$Di == 1 ) 
    
  }
  
  # add more info to dataset
  .dat$k.underlying = length(unique(.dat$study))
  .dat$k.nonaffirm.underlying = length( unique( .dat$study[ .dat$affirm == FALSE ] ) )
  
  if ( return.only.published == TRUE ) .dat = .dat[ .dat$Di == 1, ]
  
  return(.dat)
}




# always returns Cohen's d given the hack mechanisms we've chosen to implement
sim_one_study_set_stefan = function(strategy.stefan,
                                    alternative.stefan,
                                    
                                    stringent.hack,
                                    
                                    return.only.published = FALSE,
                                    is.hacked,  # separate from hack.type so that we can use same fn in each case for comparability
                                    hack.type # will be ignored if is.hacked = FALSE ("DV","optstop", "subgroup")
                                    
                                    
) {  
  
  #browser()
  
  # # test only
  # strategy.stefan = "firstsig"
  # alternative.stefan = "greater"
  # is.hacked = FALSE
  # stringent.hack = TRUE
  # return.only.published = FALSE
  # hack.type = "DV"
  
  # to force entry into while-loop below
  #sei = 100
  
  # Stefan fns routinely generate combinations of p-val and yi that imply absurdly large SEs
  #  e.g., p = 0.999 but yi = -1.04, which implies sei = 10,000
  # hackily prevent this:
  #bm
  #while ( sei > 2 | sei < 0.02 ) {
   
    
    # even if is.hacked == FALSE, we call the appropriate hack fn for the entire
    #   meta-analysis for comparability between hacked and unhacked studies, then later
    #  use the original (unhacked) stats if is.hacked = FALSE
    if (hack.type == "DV") {
      d = as.data.frame( sim.multDVhack(nobs.group = 30,
                                        nvar = 5, # default is 5
                                        r = 0.3,
                                        iter = 1,
                                        strategy = strategy.stefan,
                                        alternative = alternative.stefan,
                                        alpha = 0.05) )
    }
    
    if (hack.type == "optstop") {
      # IMPORTANT: this doesn't have strategy arg because it's inherently firstsig:
      # The dataset is evaluated row-by-row, starting with a minimum sample size of n.min. At each step, a number of observations is added to the sample, defined by the argument step and the t-test is computed. This continues until the maximum sample size specified in n.max is reached. The p-hacked p-value is defined as the first p-value that is smaller than the defined alpha level.
      
      if ( strategy.stefan != "firstsig" ) stop("hack.type optstop requires strategy = firstsig")
      
      d = as.data.frame( sim.optstop(n.min = 10,
                                     #n.max = 20,  # default
                                     n.max = 50,
                                     step = 2,
                                     #strategy = strategy.stefan,
                                     alternative = alternative.stefan,
                                     iter = 1,
                                     alpha = 0.05) )
    }
    
    
    if (hack.type == "subgroup") {
      d = as.data.frame( sim.subgroupHack(nobs.group = 30,
                                          nsubvars = 5, # default: 3
                                          strategy = strategy.stefan,
                                          alternative = alternative.stefan,
                                          alpha = 0.05,
                                          iter = 1) )
    }
    

    
    ### Post-processing that doesn't depend on hacking type ###
  
    # save the "original" (ideal draw) estimate for use with gold-std method
    d$yio = d$ds.orig
    d$pvalo = d$ps.orig
    d$seio = calc_sei(yi = d$yio, pval = d$pvalo, alternative.stefan = alternative.stefan)
  
    # choose appropriate stats as yi, sei, vi depending on whether this study is hacked
    # and get dataset into same format as in my own sim_one_study_set
    if ( is.hacked == TRUE ) {
      d = d %>% rename(pval = ps.hack,
                       yi = ds.hack) %>%
        select( !c(ps.orig, r2s.hack, r2s.orig, ds.orig) )
    } 
    if ( is.hacked == FALSE ) {
      d = d %>% rename(pval = ps.orig,
                       yi = ds.orig) %>%
        select( !c(ps.hack, r2s.hack, r2s.orig, ds.hack) )
      
    }
    
    
    # calculate sei, vi
    d$sei = calc_sei(yi = d$yi, pval = d$pval, alternative.stefan = alternative.stefan)
    sei = d$sei  # for the while-loop condition
    
  #}  # end "while ( sei > 2 | sei < 0.02 )"
  
 
  d$vi = d$sei^2
  d$tcrit = qnorm(0.975)  # here using z-approx since that's what we use to calculate sei
  
  # convenience indicators for significance and affirmative status
  d$signif = d$pval < 0.05
  d$affirm = (d$pval < 0.05 & d$yi > 0)
  
  d$is.hacked = is.hacked
  
  
  # ~~ Decide which draw to favor & publish ----
  # Stefan fn already returns only the favored draws
  # i.e., d is just a single row for this study
  #  however, need to decide whether to "report" a nonaffirm result
  if ( is.hacked == FALSE ) d$Di = 1  # because SAS only affects hacked studies
  if ( is.hacked == TRUE & stringent.hack == TRUE ) d$Di = ifelse(d$affirm == TRUE, 1, 0)
  if ( is.hacked == TRUE & stringent.hack == FALSE ) d$Di = 1
  
  if ( return.only.published == TRUE ) d = d[ d$Di == 1, ]
  
  return(d)
}





calc_sei = function(yi, pval, alternative.stefan) {
  if ( alternative.stefan == "two.sided" ) {
    # because: pval = 2 * ( 1 - pnorm( abs(yi)/sei ) )
    abs(yi) / qnorm(1 - pval/2) 
  } else if ( alternative.stefan == "greater" ) {
    # because: pval = 1 - pnorm( yi/sei )
    yi / qnorm(1 - pval) 
  } else {
    stop("calc_sei not implemented for that choice of alternative.stefan")
  }
}


# # sanity checks
# yi = -0.5
# pval = 0.03
# ( sei = calc_sei(yi = yi, pval = pval, alternative.stefan = "two.sided") )
# yi/sei  # Z-score
# expect_equal( 2 * ( 1 - pnorm( abs(yi)/sei ) ),
#               pval )
# 
# yi = -0.5
# pval = 0.8
# ( sei = calc_sei(yi = yi, pval = pval, alternative.stefan = "greater") )
# yi/sei  # Z-score
# expect_equal( 1 - pnorm( yi/sei ),
#               pval )



# DATA WRANGLING ---------------------------------------------------------------

# corrObject: something returned by correct_dataset_phack
# looks for (or makes) global object, "res"
add_method_result_row = function(repRes = NA,
                                 corrObject,
                                 methName) {
  
  # newRow = bind_cols( corrObject$metaCorr,
  #                 corrObject$sanityChecks )
  #TEMP: DON'T KEEP THE SANITY CHECKS BECAUSE CORRECT_META_PHACK2 doesn't have it
  newRow = corrObject$metaCorr
  
  newRow = newRow %>% add_column(.before = 1,
                                 methName = methName )
  
  
  # "if" condition is hacky way to deal with repRes = NA case
  if ( is.null( nrow(repRes) ) ) repRes = newRow else repRes = bind_rows(repRes, newRow)
  return(repRes)
}



# quickly look at results when running doParallel locally
srr = function(rep.res) {
  
  if( "optimx.Mhat.winner" %in% names(rep.res) ) {
    cat("\n")
    print( rep.res %>% select(method, Mhat, MLo, MHi,
                              Shat,
                              optim.converged,
                              optimx.Mhat.winner,
                              optimx.Nconvergers,
                              optimx.Pagree.of.convergers.Mhat.winner) %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  } else {
    cat("\n")
    print( rep.res %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  }
}

# SMALL GENERIC HELPERS ---------------------

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}


# check CI coverage
covers = function( truth, lo, hi ) {
  return( (lo <= truth) & (hi >= truth) )
}

# get names of dataframe containing a string
namesWith = function(pattern, dat){
  names(dat)[ grepl(pattern = pattern, x = names(dat) ) ]
}


# quick length(unique)
nuni = function(x) {
  length(unique(x))
}

# (re-)install package AND its dependencies
# useful for stupid rstan issues in which rstan itself it UTD but not its dependencies
# https://stackoverflow.com/questions/21010705/update-a-specific-r-package-and-its-dependencies
instPkgPlusDeps <- function(pkg, install = FALSE,
                            which = c("Depends", "Imports", "LinkingTo"),
                            inc.pkg = TRUE) {
  stopifnot(require("tools")) ## load tools
  ap <- available.packages() ## takes a minute on first use
  ## get dependencies for pkg recursively through all dependencies
  deps <- package_dependencies(pkg, db = ap, which = which, recursive = TRUE)
  ## the next line can generate warnings; I think these are harmless
  ## returns the Priority field. `NA` indicates not Base or Recommended
  pri <- sapply(deps[[1]], packageDescription, fields = "Priority")
  ## filter out Base & Recommended pkgs - we want the `NA` entries
  deps <- deps[[1]][is.na(pri)]
  ## install pkg too?
  if (inc.pkg) {
    deps = c(pkg, deps)
  }
  ## are we installing?
  if (install) {
    install.packages(deps)
  }
  deps ## return dependencies
}

# example
# instPkgPlusDeps("fields")


# CLUSTER FNS ---------------------------------------------------------------

# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
sbatch_skeleton <- function() {
  return(
    "#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=JOBNAME
#################  
#a file for job output, you can check job progress
#SBATCH --output=OUTFILE
#################
# a file for errors from the job
#SBATCH --error=ERRORFILE
#################
#time you think you need; default is one hour
#SBATCH --time=JOBTIME
#################
#quality of service; think of it as job priority
#SBATCH --qos=QUALITY
#################
#submit to both owners and normal partition
#SBATCH -p normal,owners,qsu
#################
#number of nodes you are requesting
#SBATCH --nodes=NODENUMBER
#################
#memory per node; default is 4000 MB
#SBATCH --mem=MEMPERNODE
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=MAILTYPE
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=USER_EMAIL
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=TASKS_PER_NODE
#################
#SBATCH --cpus-per-task=CPUS_PER_TASK
#now run normal batch commands

ml load v8
ml load R/4.2.0
ml load jags
R -f PATH_TO_R_SCRIPT ARGS_TO_R_SCRIPT")
}



generateSbatch <- function(sbatch_params,
                           runfile_path = NA,
                           run_now = F) {
  
  #sbatch_params is a data frame with the following columns
  #jobname: string, specifies name associated with job in SLURM queue
  #outfile: string, specifies the name of the output file generated by job
  #errorfile: string, specifies the name of the error file generated by job
  #jobtime: string in hh:mm:ss format, max (maybe soft) is 48:00:00 
  #specifies the amoung of time job resources should be allocated
  #jobs still running after this amount of time will be aborted
  #quality: kind of like priority, normal works
  #node_number, integer: the number of nodes (computers w/16 cpus each) to allocate 
  #mem_per_node, integer: RAM, in MB, to allocate to each node
  #mailtype, string: ALL, BEGIN, END, FAIL: what types of events should you be notified about via email
  #user_email string: email address: email address to send notifications
  #tasks_per_node: integer, number of tasks, you should probably use 1
  #cpus_per_task: integer, 1-16, number of cpus to use, corresponds to number of available cores per task
  #path_to_r_script: path to r script on sherlock
  #args_to_r_script: arguments to pass to r script on command line
  #write_path: where to write the sbatch file
  #server_sbatch_path: where sbatch files will be stored on sherlock
  #runfile_path is a string containing a path at which to write an R script that can be used to run
  #the batch files generated by this function. 
  #if NA, no runfile will be written
  #run_now is a boolean specifying whether batch files should be run as they are generated
  
  sbatches <- list()
  if (!is.na(runfile_path)) {
    outfile_lines <- c(paste0("# Generated on ",  Sys.time()))
  }
  for (sbatch in 1:nrow(sbatch_params) ) {
    gen_batch <- sbatch_skeleton()
    #set job name
    if (is.null(sbatch_params$jobname[sbatch])) { 
      gen_batch <- gsub("JOBNAME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBNAME", sbatch_params$jobname[sbatch], gen_batch) 
    }
    #set outfile name
    if (is.null(sbatch_params$outfile[sbatch])) { 
      gen_batch <- gsub("OUTFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("OUTFILE", sbatch_params$outfile[sbatch], gen_batch) 
    }
    #set errorfile name
    if (is.null(sbatch_params$errorfile[sbatch])) { 
      gen_batch <- gsub("ERRORFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ERRORFILE", sbatch_params$errorfile[sbatch], gen_batch) 
    }
    #set jobtime
    if (is.null(sbatch_params$jobtime[sbatch])) { 
      gen_batch <- gsub("JOBTIME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBTIME", sbatch_params$jobtime[sbatch], gen_batch) 
    }
    #set quality
    if (is.null(sbatch_params$quality[sbatch])) { 
      gen_batch <- gsub("QUALITY", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("QUALITY", sbatch_params$quality[sbatch], gen_batch) 
    }
    #set number of nodes
    if (is.null(sbatch_params$node_number[sbatch])) { 
      gen_batch <- gsub("NODENUMBER", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("NODENUMBER", sbatch_params$node_number[sbatch], gen_batch) 
    }
    #set memory per node
    if (is.null(sbatch_params$mem_per_node[sbatch])) { 
      gen_batch <- gsub("MEMPERNODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MEMPERNODE", sbatch_params$mem_per_node[sbatch], gen_batch) 
    }
    #set requested mail message types
    if (is.null(sbatch_params$mailtype[sbatch])) { 
      gen_batch <- gsub("MAILTYPE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MAILTYPE", sbatch_params$mailtype[sbatch], gen_batch) 
    }
    #set email at which to receive messages
    if (is.null(sbatch_params$user_email[sbatch])) { 
      gen_batch <- gsub("USER_EMAIL", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("USER_EMAIL", sbatch_params$user_email[sbatch], gen_batch) 
    }
    #set tasks per node
    if (is.null(sbatch_params$tasks_per_node[sbatch])) { 
      gen_batch <- gsub("TASKS_PER_NODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("TASKS_PER_NODE", sbatch_params$tasks_per_node[sbatch], gen_batch) 
    }
    #set cpus per task
    if (is.null(sbatch_params$cpus_per_task[sbatch])) { 
      gen_batch <- gsub("CPUS_PER_TASK", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("CPUS_PER_TASK", sbatch_params$cpus_per_task[sbatch], gen_batch) 
    }
    #set path to r script
    if (is.null(sbatch_params$path_to_r_script[sbatch])) { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", sbatch_params$path_to_r_script[sbatch], gen_batch) 
    }
    #set args to r script
    if (is.null(sbatch_params$args_to_r_script[sbatch])) { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", sbatch_params$args_to_r_script[sbatch], gen_batch) 
    }
    
    #write batch file
    if (is.null(sbatch_params$write_path[sbatch])) { 
      cat(gen_batch, file = paste0("~/sbatch_generated_at_", gsub(" |:|-", "_", Sys.time()) ), append = F)
    } else { 
      cat(gen_batch, file = sbatch_params$write_path[sbatch], append = F)
    }
    
    if (!is.na(sbatch_params$server_sbatch_path[sbatch])) {
      outfile_lines <- c(outfile_lines, paste0("system(\"sbatch ", sbatch_params$server_sbatch_path[sbatch], "\")"))
    } 
    sbatches[[sbatch]] <- gen_batch
  }
  if (!is.na(runfile_path)) {
    cat(paste0(outfile_lines, collapse = "\n"), file = runfile_path)
  }
  if(run_now) { system(paste0("R -f ", runfile_path)) } 
  
  return(sbatches)
}


# looks at results files to identify sbatches that didn't write a file
# .max.sbatch.num: If not passed, defaults to largest number in actually run jobs.

sbatch_not_run = function(.results.singles.path,
                          .results.write.path,
                          .name.prefix,
                          .max.sbatch.num = NA ) {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # extract job numbers
  sbatch.nums = as.numeric( unlist( lapply( strsplit( keepers, split = "_"), FUN = function(x) x[5] ) ) )
  
  # check for missed jobs before the max one
  if ( is.na(.max.sbatch.num) ) .max.sbatch.num = max(sbatch.nums)
  all.nums = 1 : .max.sbatch.num
  missed.nums = all.nums[ !all.nums %in% sbatch.nums ]
  
  # give info
  print( paste("The max job number is: ", max(sbatch.nums) ) )
  print( paste( "Number of jobs that weren't run: ",
                ifelse( length(missed.nums) > 0, length(missed.nums), "none" ) ) )
  
  if( length(missed.nums) > 0 ) {
    setwd(.results.write.path)
    write.csv(missed.nums, "missed_job_nums.csv")
  }
  
  return(missed.nums)
  
}

# FN: STITCH RESULTS FILES -------------------------------------

# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path, .results.stitched.write.path=.results.singles.path,
                        .name.prefix, .stitch.file.name="stitched_model_fit_results.csv") {
  
  # .results.singles.path = "/home/groups/manishad/MRM/sim_results/long"
  # .results.stitched.write.path = "/home/groups/manishad/MRM/sim_results/overall_stitched"
  # .name.prefix = "long_results"
  # .stitch.file.name="stitched.csv"
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # grab variable names from first file
  names = names( read.csv(keepers[1] )[-1] )
  
  # read in and rbind the keepers
  tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )
  s <- do.call(rbind, tables)
  
  names(s) = names( read.csv(keepers[1], header= TRUE) )
  
  if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
  write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
  return(s)
}


# quickly look at results from job #1

res1 = function() {
  setwd("/home/groups/manishad/SAPH/long_results")
  rep.res = fread("long_results_job_1_.csv")
  srr()
  
  cat("\nErrors by method:" )
  print( rep.res %>% group_by(method) %>%
           summarise(prop.error = mean( overall.error != "" ) ) )
  
  #table(rep.res$overall.error)
  
  cat("\n\nDim:", dim(rep.res))
  cat("\n\nReps completed:", nrow(rep.res)/nuni(rep.res$method))
}




# TRASH -------------------------------------

# # 2022-3-19: I NO LONGER THINK IT MAKES SENSE TO USE THIS PLOT WITH VARIABLE SEI'S. 
# #  BUT SAVE IT B/C USEFUL WHEN SEI'S ARE THE SAME.
# #  Since the sei's differ, so do the cutpoints, and even if we plot the full normal dist, 
# #  its height doesn't align properly with the density when the density is truncated,
# #  giving the impression of a poor fit.
# # Args:
# #  - d: dataset with var names "yi", "vi", "affirm"
# #  - Mhat: RTMA estimate of underlying mean effect size
# #  - Shat: Same for heterogeneity
# #  - showAffirms: should it show all studies, even affirms?
# 
# # Returned plot:
# #  - black line = LOESS density of nonaffirms
# #  - red line = MLE from RTMA (parametric counterpart to the above)
# #  - blue line = LOESS density of all tstats (including affirms)
# 
# # IMPORTANT:
# # Note that truncation point shown in plots is only approximate because it’s set to 1.96 for all Zi.tilde, whereas actually each Zi.tilde has its own trunc point (see Zi_tilde_cdf).
# plot_trunc_densities_RTMA = function(d,
#                                      Mhat,
#                                      Shat,
#                                      showAffirms = FALSE) {
#   
#   # #TEST ONLY
#   # d = dp
#   # showAffirms = FALSE
#   # Mhat = 0.45  # FAKE for now
#   # Shat = 0.10
#   
#   # add Z-scores
#   # these are standardized ACROSS studies using the ESTIMATED mean and heterogeneity
#   # so they should look truncated N(0,1)
#   d$Zi.tilde = (d$yi - Mhat) / sqrt(Shat^2 + d$vi)
#   
#   # already has affirmative indicator
#   dn = d %>% filter(affirm == FALSE)
#   
#   
#   xmin = floor(min(dn$Zi.tilde))
#   xmax = ceiling(max(dn$Zi.tilde))
#   
#   p = ggplot(data = data.frame(x = c(xmin, 3)),
#              aes(x)) +
#     
#     geom_vline(xintercept = 0,
#                lwd = 1,
#                color = "gray") +
#     
#     # estimated density of estimates
#     geom_density( data = dn,
#                   aes(x = Zi.tilde),
#                   adjust = .3 ) +
#     
#     # estimated density from meta-analysis
#     # stat_function( fun = dnorm,
#     #                n = 101,
#     #                args = list( 
#     #                             mean = 0,
#     #                             sd = 1 ),
#     #                #aes(y = .25 * ..count..),  # doesn't work
#     #                lwd = 1.2,
#     #                color = "red") +
#     stat_function( fun = dtrunc,
#                    n = 101,
#                    args = list( spec = "norm",
#                                 mean = 0,
#                                 sd = 1,
#                                 b = qnorm(.975) ),
#                    #aes(y = .25 * ..count..),  # doesn't work
#                    lwd = 1.2,
#                    color = "red") +
#     #     
#     
#     ylab("") +
#     #scale_x_continuous( breaks = seq(xmin, 3, 0.5)) +
#     xlab("Across-study Z-score using (Mhat, Shat)") +
#     theme_minimal() +
#     scale_y_continuous(breaks = NULL) +
#     theme(text = element_text(size=16),
#           axis.text.x = element_text(size=16))
#   
#   
#   # also show density of all t-stats, not just the nonaffirms
#   if ( showAffirms == TRUE ) {
#     p = p + geom_density( data = d,
#                           aes(x = Zi.tilde),
#                           color = "blue",
#                           adjust = .3 )
#   }
#   
#   
#   return(p)
#   
# }
# # OLDER VERSION - MAYBE SAVE?
# # plot empirical data
# 
# # .obj: object returned by correct_meta_phack2
# # showAffirms: should it show all studies, even affirms?
# # black line = LOESS density of nonaffirms
# # red line = MLE from RTMA (parametric counterpart to the above)
# # blue line = LOESS density of all tstats (including affirms)
# plot_trunc_densities = function(.obj,
#                                 showAffirms = FALSE) {
#   
#   # already has affirmative indicator
#   d = .obj$data
#   dn = d[d$affirm == FALSE,]
#   
#   tstatMeanMLE = .obj$sanityChecks$tstatMeanMLE
#   tstatVarMLE = .obj$sanityChecks$tstatVarMLE
#   
#   xmin = floor(min(dn$tstat))
#   xmax = ceiling(max(dn$tstat))
#   
#   p = ggplot(data = data.frame(x = c(xmin, 3)),
#              aes(x)) +
#     
#     geom_vline(xintercept = 0,
#                lwd = 1,
#                color = "gray") +
#     
#     # geom_vline(xintercept = tstatMeanMLE,
#     #            lty = 2,
#     #            lwd = 1,
#     #            color = "red") +
#     
#     
#     # estimated density of estimates
#     geom_density( data = dn,
#                   aes(x = tstat),
#                   adjust = .3 ) + 
#     
#     
#     
#     # estimated density from meta-analysis
#     stat_function( fun = dtrunc,
#                    n = 101,
#                    args = list( spec = "norm",
#                                 mean = tstatMeanMLE,
#                                 sd = sqrt(tstatVarMLE),
#                                 b = .obj$crit),
#                    #aes(y = .25 * ..count..),  # doesn't work
#                    lwd = 1.2,
#                    color = "red") +
#     
#     
#     
#     ylab("") +
#     #scale_x_continuous( breaks = seq(xmin, 3, 0.5)) +
#     xlab("t-stat") +
#     theme_minimal() +
#     scale_y_continuous(breaks = NULL) +
#     theme(text = element_text(size=16),
#           axis.text.x = element_text(size=16))
#   
#   
#   # also show density of all t-stats, not just the nonaffirms
#   if ( showAffirms == TRUE ) {
#     p = p + geom_density( data = d,
#                           aes(x = tstat),
#                           color = "blue",
#                           adjust = .3 )
#   }
#   
#   
#   return(p)
#   
# }

# OLD VERSION (for sanity-checking the one below)

# 2022-3-12
# nonaffirms only
# Zi_tilde_cdf_OLD = function(x, .SE, .Shat) {
#   
#   # calculate cutpoint for EACH Zi.tilde
#   # **reasoning:
#   #  we observe a truncated sample st yi > 1.96*SE
#   # therefore Zi.tilde = (yi - mu) / ( sqrt(Shat^2 + SE^2) ) > (1.96*SE - mu) / ( sqrt(Shat^2 + SE^2) )
#   # and we know that Zi.tilde ~ N(0,1) prior to truncation
#   Zi.tilde.crit = ( qnorm(.975) * .SE ) / sqrt(.Shat^2 + .SE^2)
#   
#   ptruncnorm(q = x,
#              a = -Inf,
#              b = Zi.tilde.crit,
#              mean = 0,
#              sd = 1)
# }

# # 2022-3-12
# # fit diagnostics
# # get CDF of (non-iid) marginal Z-scores (Zi.tilde)
# #  given a fitted Shat
# # .affirm: VECTOR with same length as x for affirm status
# #  including the affirms is useful for 2PSM
# Zi_tilde_cdf = function(.Zi.tilde,
#                         .SE,
#                         .Shat,
#                         .affirm) {
#   
#   
#   #if ( length(.Zi.tilde) > 1 ) stop("Length of .Zi.tilde must be 1")
#   if ( length(.affirm) != length(.Zi.tilde) ) stop(".affirm must have same length as x")
#   
#   # calculate cutpoint for this Zi.tilde
#   # **reasoning:
#   #  we observe a truncated sample st yi > 1.96*SE
#   # therefore Zi.tilde = (yi - mu) / ( sqrt(Shat^2 + .SE^2) ) > (1.96*SE - mu) / ( sqrt(Shat^2 + .SE^2) )
#   # and we know that Zi.tilde ~ N(0,1) prior to truncation
#   Zi.tilde.crit = ( qnorm(.975) * .SE ) / sqrt(.Shat^2 + .SE^2)
#   
#   
#   dat = data.frame(Zi.Tilde = .Zi.tilde,
#                    Zi.Tilde.Crit = Zi.tilde.crit,
#                    Affirm = .affirm)
#   
#   # if ( .affirm == FALSE ) expect_equal( .Zi.tilde < Zi.tilde.crit, TRUE )
#   # if ( .affirm == FALSE ) expect_equal( .Zi.tilde < Zi.tilde.crit, TRUE )
#   
#   # if ( Affirm == FALSE ) {
#   #   return( ptruncnorm(q = Zi.Tilde,
#   #                      a = -Inf,
#   #                      b = Zi.Tilde.Crit,
#   #                      mean = 0,
#   #                      sd = 1) )
#   # } else if ( Affirm == TRUE ) {
#   #   return( ptruncnorm(q = Zi.Tilde,
#   #                      a = Zi.Tilde.Crit,
#   #                      b = Inf,
#   #                      mean = 0,
#   #                      sd = 1) )
#   # }
#   
#   dat$cdfi = NA
#   
#   if ( any(dat$Affirm == FALSE) ) {
#     dat$cdfi[ dat$Affirm == FALSE ] = ptruncnorm(q = dat$Zi.Tilde[ dat$Affirm == FALSE ],
#                                                  a = -Inf,
#                                                  b = dat$Zi.Tilde.Crit[ dat$Affirm == FALSE ],
#                                                  mean = 0,
#                                                  sd = 1)
#   }
#   
#   if ( any(dat$Affirm == TRUE) ) {
#     dat$cdfi[ dat$Affirm == TRUE ] = ptruncnorm(q = dat$Zi.Tilde[ dat$Affirm == TRUE ],
#                                                 a = dat$Zi.Tilde.Crit[ dat$Affirm == TRUE ],
#                                                 b = Inf,
#                                                 mean = 0,
#                                                 sd = 1)
#   }
#   
#   return(dat$cdfi)
#   
#   
#   # 
#   # dat = dat %>% rowwise() %>%
#   #   mutate( cdfi = function(Zi.Tilde, Affirm){
#   #     if ( Affirm == FALSE ) {
#   #       return( ptruncnorm(q = Zi.Tilde,
#   #                          a = -Inf,
#   #                          b = Zi.Tilde.Crit,
#   #                          mean = 0,
#   #                          sd = 1) )
#   #     } else if ( Affirm == TRUE ) {
#   #       return( ptruncnorm(q = Zi.Tilde,
#   #                          a = Zi.Tilde.Crit,
#   #                          b = Inf,
#   #                          mean = 0,
#   #                          sd = 1) )
#   #     }
#   #     
#   #     
#   #   } )
#   
# }


# ### Sanity checks ###
# 
# # sanity-check nonaffirmatives (taken from Kvarven's Belle meta-analysis)
# Zi.tilde = c(-1.19493855966001, -0.330782431293096, 0.18135714022493, -0.355495378590117)
# sei = c(0.183, 0.169, 0.269999994444444, 0.222)
# Shat = 0.23  # from 2PSM
# cdfi.mine = ptruncnorm(q = Zi.tilde,
#                            a = -Inf,
#                            b = ( qnorm(.975) * sei ) / sqrt(Shat^2 + sei^2),
#                            mean = 0,
#                            sd = 1)
# 
# 
# cdfi = Zi_tilde_cdf(.Zi.tilde = Zi.tilde,
#                     .SE = sei,
#                     .Shat = Shat,
#                     .affirm = rep(FALSE, length(Zi.tilde)))
# expect_equal( cdfi.mine, cdfi)
# 
# # sanity-check affirmatives (also from Belle)
# Zi.tilde = c(2.10452951219512, 2.28191480814403, 4.29142844422158, 2.33576635036496, 
#              2.707112988856, 2.86885247443619, 3.10843366971151, 2.37931028904956, 
#              2.6049383759669, 2.06417123549919, 2.031579)
# sei = c(0.287, 0.188000002659574, 0.175000002857143, 0.137, 0.23900000209205, 
#         0.304999998360656, 0.249000002008032, 0.261000001915709, 0.242999997942387, 
#         0.186999994652406, 0.19)
# Shat = 0.23
# cdfi.mine = ptruncnorm(q = Zi.tilde,
#                        a = ( qnorm(.975) * sei ) / sqrt(Shat^2 + sei^2),
#                        b = Inf,
#                        mean = 0,
#                        sd = 1)
# 
# 
# cdfi = Zi_tilde_cdf(.Zi.tilde = Zi.tilde,
#                     .SE = sei,
#                     .Shat = Shat,
#                     .affirm = rep(TRUE, length(Zi.tilde)))
# expect_equal( cdfi.mine, cdfi)



# # 2022-3-12
# # test for RTMA fit
# # yi, sei: can be for all published studies or for just nonaffirms or just affirms
# #  including the affirms is useful for 2PSM but not RTMA
# 
# my_ks_test_RTMA = function(yi,
#                            sei,
#                            Mhat,
#                            Shat) {
#   
#   if ( is.na(Mhat) | is.na(Shat) ) return(NA)
#   
#   affirm = (yi/sei) > qnorm(.975)
#   
#   # retain only nonaffirmatives
#   # nonaffirm = yi/sei < qnorm(.975)
#   # yi = yi[ nonaffirm == TRUE ]
#   # sei = sei[ nonaffirm == TRUE ]
#   
#   Zi.tilde = (yi - Mhat) / sqrt(Shat^2 + sei^2)
#   
#   res = ks.test(Zi.tilde, function(x) Zi_tilde_cdf(.Zi.tilde = Zi.tilde,
#                                                    .SE = sei,
#                                                    .Shat = Shat,
#                                                    .affirm = affirm) )
#   res$p.value
#   
# }




# # 2022-3-12
# # test for RTMA fit
# # yi, sei: can be for all published studies or for just nonaffirms; 
# #  fn will automatically retain only nonaffirms
# my_ks_test_RTMA = function(yi,
#                            sei,
#                            Mhat,
#                            Shat) {
#   
#   # retain only nonaffirmatives
#   nonaffirm = yi/sei < qnorm(.975)
#   yi = yi[ nonaffirm == TRUE ]
#   sei = sei[ nonaffirm == TRUE ]
#   
#   Zi.tilde = (yi - Mhat) / sqrt(Shat^2 + sei^2)
#   
#   res = ks.test(Zi.tilde, function(x) Zi_tilde_cdf(x,
#                                                    .SE = sei,
#                                                    .Shat = Shat) )
#   res$p.value
#   
# }
