
# NOTES ---------------------------------------------------------------

# "draw" always refers to within-study draw
# "study set" refers to all draws, observed and observed, for a given study


# ESTIMATION METHOD FNS ----------------------------------------------

# ~ Estimation Methods Structured for Use Inside run_method_safe --------

# Because these fns are run inside run_method_safe, the latter will handle editing rep.res
#  All of these fns should take get.CIs as an argument and return CIs as c(NA, NA) if not wanted
# Fns in this category need to return a dataframe with the below structure, although it's okay if they don't return all of these names since run_method_safe will handle that. Note that this is a LIST containing a dataframe called "stats", not just the dataframe; this allows easy extension in case you want to return other objects, like model objects.

# .yi: published point estimates
# .sei: their SEs
# .tcrit: critical values on t or z scale for each study; can just use qnorm(.975) by default
# .Mu.start: optimizer starting value for meta-analysis Mu
# .Tt.start: optimizer starting value for meta-analysis tau
# .stan.adapt_delta: passed to rstan
# .stan.maxtreedepth: same
#  we should later use ellipsis to allow passing arbitrary args to rstan
estimate_jeffreys = function(.yi,
                             .sei,
                             
                             .Mu.start,
                             .Tt.start,
                             .stan.adapt_delta = 0.8,
                             .stan.maxtreedepth = 10 ) {
  
  # stan.model (used later) is compiled OUTSIDE this fn in doParallel to avoid 
  #  issues with nodes competing with one another
  
  
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
    
    cat( paste("\n estimate_jeffreys flag 2: about to call sampling") )
    
    post = sampling(stan.model,
                    cores = 1,
                    refresh = 0,
                    data = list( k = length(.yi),
                                 sei = .sei,
                                 y = .yi ),
                    
                    #iter = p$stan.iter,   
                    control = list(max_treedepth = .stan.maxtreedepth,
                                   adapt_delta = .stan.adapt_delta),
                    
                    init = init.fcn)
    
    
  }, warning = function(condition){
    stan.warned <<- 1
    stan.warning <<- condition$message
  } )
  
  cat( paste("\n estimate_jeffreys flag 3: about to call postSumm") )
  postSumm = summary(post)$summary
  if (is.null(postSumm)) stop("In stan, postSumm is null")
  
  # pull out best iterate to pass to MAP optimization later
  ext = rstan::extract(post) # a vector of all post-WU iterates across all chains
  best.ind = which.max(ext$log_post)  # single iterate with best log-posterior should be very close to MAP
  
  
  # posterior means, posterior medians, modes, and max-LP iterate
  Mhat = c( postSumm["mu", "mean"],
            median( rstan::extract(post, "mu")[[1]] ),
            ext$mu[best.ind] )
  
  Shat = c( postSumm["tau", "mean"],
            median( rstan::extract(post, "tau")[[1]] ),
            ext$tau[best.ind] )
  
  # sanity check
  #expect_equal( Mhat[1], mean( rstan::extract(post, "mu")[[1]] ) )
  
  
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
    
    # this will use same CI limits for all pt estimates
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



# FNS FOR FUTURE R PACKAGE ----------------------------

# structured as in phacking pkg

get_lprior <- function(mu, tau, sei) {
  e_fisher_i <- function(se) {
    si <- sqrt(tau ^ 2 + se ^ 2)
    
    kmm <- -si ^ (-2)
    kms <- 0
    kss <- -2 * tau ^ 2 * si ^ (-4)
    
    matrix(c(-kmm, -kms, -kms, -kss), nrow = 2, ncol = 2)
  }
  
  e_fisher <- purrr::map(sei, e_fisher_i) |> purrr::reduce(`+`)
  log(sqrt(det(e_fisher)))
}

get_nll <- function(mu, tau, yi, sei) {
  si <- sqrt(tau ^ 2 + sei ^ 2)
  sum(log(si * sqrt(2 * pi)) + 0.5 * si ^ (-2) * (yi - mu) ^ 2)
}

nlpost <- function(mu, tau, yi, sei) {
  joint_nll <- get_nll(mu, tau, yi, sei) # negative log-likelihood
  joint_lprior <- get_lprior(mu, tau, sei) # log-prior
  joint_nll - joint_lprior # log-posterior
}

mle_params <- function(mu_start, tau_start, yi, sei) {
  nlpost_fun <- function(mu, tau) nlpost(mu, tau, yi, sei)
  stats4::mle(minuslogl = nlpost_fun,
              start = list(mu = mu_start, tau = tau_start),
              method = "Nelder-Mead")
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



# DATA SIMULATION ---------------------------------------------------------------


# - Mu: overall mean for meta-analysis (SMD if Ytype = "cont-SMD"; log-RR if Ytype = "bin-RR"; log-OR if Ytype = "log-OR")
# - t2a: across-study heterogeneity (NOT total heterogeneity)
# - true.dist: "norm" or "expo"
# - muN, minN: mean and lower limit of uniform dist from which to draw sample sizes
# - Ytype: "cont-SMD" (SMD effect size), "bin-RR" (RR effect size), "bin-OR" (OR effect size)
# - p0: P(Y=0 | X=0); only needed if Ytype is binary 
sim_meta = function(k.pub,
                    Mu,  
                    t2a,  
                    true.dist,
                    
                    # within-study parameters
                    N.expr,
                    Ytype,
                    p0) {
  
  
  # collect arguments
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k.pub")]
  
  #browser()

  for( i in 1:k.pub ) {
    
    newRow = do.call( sim_one_study, .args )
    
    # add study ID
    newRow = newRow %>% add_column( .before = 1,
                                      study = i )
    
    if ( i == 1 ) .dat = newRow else .dat = rbind( .dat, newRow )
  }
  
  # add more info to dataset
  .dat$mean_pY = meanNA(.dat$pY)
  .dat$mean_yi = meanNA(.dat$yi)
  
  return(.dat)
}

# test
if (FALSE) {
  # example: continuous Y
  d = sim_meta(k.pub = 10,
               
               Mu = 0.5,  
               t2a = 0.1^2,  
               true.dist = "norm",
               
               # within-study parameters
               # muN = 1000,
               # minN = 1000,
               N.expr = "round( runif(n = 1, min = 40, max = 400) )",
               Ytype = "cont-SMD",
               p0 = NA)
  
  # example: binary Y
  # evil scen 105
  d = sim_meta(  k.pub = 100,
                 t2a = 0.0001,
                 Mu = 0,
                 true.dist = "norm",
                 p0 = 0.05,
                 Ytype = "bin-OR",
                 N.expr = "40" )
  hist(d$sei)
  
  mean(d$pY)
  mean(d$pY0) # should match p0
}



# ~ Simulate a single study ----------------- 

# see sim_meta for args
sim_one_study = function( Mu,  # overall mean for meta-analysis
                          t2a,  # across-study heterogeneity
                          true.dist,
                          
                          # within-study sample size parameters
                          N.expr, 
                          
                          sd.w = 1,  # within-study SD(Y|X); only needed for cont outcome
                          
                          Ytype,  
                          p0 = NULL # P(Y | X=0); only needed for binary outcome
) {  
  
  # for testing
  if (FALSE){
    true.dist = "expo"
    Mu = 0.5
    t2a = 0
    muN = minN = 40
    Ytype = "bin-OR"
    sd.w = 1
    p0 = 0.01
  }
  
  # ~~ Mean for this study set -------------------------------------------------
  
  if( !true.dist %in% c("norm", "expo") ) stop("true.dist not recognized")
  
  if ( true.dist == "norm" ){
    mui = Mu + rnorm(mean = 0,
                     sd = sqrt(t2a),
                     n = 1)
  }
  
  if ( true.dist == "expo" ){
    # set the rate so the heterogeneity is correct
    mui = rexp( n = 1, rate = sqrt(1/t2a) )
    # now the mean is sqrt(t2a) rather than Mu
    # shift to have the correct mean (in expectation)
    mui = mui + ( Mu - sqrt(t2a))
  }
  
  # simulate total N for this study
  #if ( muN < minN ) stop("Should not have muN < minN")
  #N = round( runif( n = 1, min = minN, max = minN + 2*( muN - minN ) ) ) # draw from uniform centered on muN
  N = eval( parse( text = N.expr ) )
  muN = N  #@TEMP FOR USING N.EXPR INSTEAD OF MUN, MINN (needed for sanchecks)
  
  # ~~ Simulate individual subject data -------------------------------------------------
  
  # as in MRM helper code
  if ( Ytype == "cont-SMD" ) {
    # group assignments
    X = c( rep( 0, N/2 ), rep( 1, N/2 ) )
    
    ### Continuous Y ###
    # 2-group study of raw mean difference with means 0 and Mi in each group
    # and same SD
    Y = c( rnorm( n = N/2, mean = 0, sd = sd.w ),
           rnorm( n = N/2, mean = mui, sd = sd.w ) )
    
    # calculate ES for this study using metafor (see Viechtbauer "Conducting...", pg 10)
    ES = escalc( measure="SMD",   
                 n1i = N/2, 
                 n2i = N/2,
                 m1i = mean( Y[X==1] ),
                 m2i = mean( Y[X==0] ),
                 sd1i = sd( Y[X==1] ),
                 sd2i = sd( Y[X==0] ) ) 
    
    # only here to be consistent with binary Y case below
    pY = pY1 = pY0 = nY1 = nY1_theory = nY0 = nY0_theory = NA
  }
  
  
  # similar to Metasens helper code
  if ( Ytype %in% c("bin-RR", "bin-OR") ) {
    
    if ( is.null(p0) ) stop("Must specify p0 for binary Y")
    
    # group assignments
    X = c( rep( 0, N/2 ), rep( 1, N/2 ) )
    
    ### Binary Y; odds ratio ###
    if (Ytype == "bin-RR") {
      
      # check that args are ok
      if ( p0 * exp(mui) > 1 ) stop("Theoretical P(Y=1 | X=1) > 1. Adjust p0 or Mu.")
      if ( p0 * exp(mui) < 0 ) stop("Theoretical P(Y=1 | X=1) < 0. Adjust p0 or Mu.")
      
      linpred = log(p0) + mui*X  # mui is already on log scale
      # exp here because log-RR model
      Y = rbinom( size=1, n=N, prob=exp(linpred) ) 
      
      # sanity check to be returned
      nY0_theory = p0 * (muN/2)
      nY1_theory = p0 * exp(mui) * (muN/2)
      
      # sanity check
      if (FALSE){
        coef( glm(Y ~ X, family=binomial(link = "log")) )[["X"]]; mui
        mean(Y[X==0]); p0
        mean(Y[X==1]); p0 * exp(mui)
      }
    
    ### Binary Y; risk ratio ###
    } else if (Ytype == "bin-OR") {
      
      # no need to check that args are ok as above, since expit in [0,1]
      
      linpred = logit(p0) + mui*X 
      Y = rbinom( size=1, n=N, prob=expit(linpred) ) 
      
      # sanity check to be returned
      nY0_theory = p0 * (muN/2)
      nY1_theory = expit( logit(p0) + mui ) * (muN/2)

      # sanity check
      if (FALSE){
        coef( glm(Y ~ X, family=binomial(link = "logit")) )[["X"]]; mui
        mean(Y[X==0]); p0
        mean(Y[X==1]); expit( logit(p0) + mui )
      }
    }
    
    
    # calculate deaths (Y=1) and sample sizes in each group
    n1 = sum(X)  # number deaths among X=1
    n0 = length(X) - sum(X)
    y1 = sum(Y[X==1])  # number of deaths in Tx group
    y0 = sum(Y[X==0])  # number of deaths in control group
    
    # calculate log-RR for this study using metafor (see Viechtbauer "Conducting...", pg 10)
    if (Ytype == "bin-RR") {
      ES = escalc( measure="RR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0 )  # returns on log scale
      
    } else if (Ytype == "bin-OR") {
      
      ES = escalc( measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0 )  # returns on log scale
    }
    
    # sanity checks: summary stats
    pY = mean(Y)
    pY1 = mean(Y[X==1])
    pY0 = mean(Y[X==0])
    nY1 = y1
    nY0 = y0
  }
  
  yi = ES$yi
  vi = ES$vi
  sei = sqrt(vi)
  
  # ~~ One row for meta-analytic dataset  -------------------------------------------------
  d = data.frame( yi = yi, 
                  vi = vi, 
                  sei = sei, 
                  N = N,
                  
                  pY = pY,
                  pY1 = pY1,
                  pY0 = pY0,
                  
                  nY1 = nY1,
                  nY1_theory = nY1_theory,
                  
                  nY0 = nY0,
                  nY0_theory = nY0_theory)
  
  return(d)
  
}



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

# take the logit of a probability, but truncate
#  to avoid infinities
truncLogit <- function(p) {
  p[p==0] = 0.001
  p[p==1] = 0.999
  log(p/(1-p))
}


logit <- function(p) {
  log(p/(1-p))
}


expit = function(x) {
  exp(x) / (1 + exp(x))
}

# calculate I^2 from t^2 and N
I2 = function(t2, N) {
  t2 / (t2 + 4/N)
}

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}

# quick mean with NAs removed
medianNA = function(x){
  median(x, na.rm = TRUE)
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
  if ( is.na(.max.sbatch.num) ) .max.sbatch.num = max(sbatch.nums, na.rm = TRUE)
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
