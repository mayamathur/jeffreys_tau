
# PRELIMINARIES -----------------------------

# rm( list = ls() )

# are we running locally?
run.local = FALSE
# run.local = TRUE

# should we set scen params interactively on cluster?
interactive.cluster.run = FALSE

# ~~ Packages -----------------------------------------------
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
           "rma.exact",
           "CompQuadForm",
           "bayesmeta",
           "phacking")  # note: to reinstall this one, need ml load jags

# to install everything
# lapply(toLoad, install.packages)

if ( run.local == TRUE | interactive.cluster.run == TRUE ) toLoad = c(toLoad, "here")



# SET UP FOR CLUSTER OR LOCAL RUN ------------------------------

# ~~ Local Run ----------------------------------------
if (run.local == FALSE) {
  
  # load command line arguments
  args = commandArgs(trailingOnly = TRUE)
  
  cat("\n\n args received from sbatch file:", args)
  
  jobname = args[1]
  scen = args[2]  # this will be a number
  
  # load packages with informative messages if one can't be installed
  # **Common reason to get the "could not library" error: You did ml load R/XXX using an old version
  any.failed = FALSE
  for (pkg in toLoad) {
    
    cat( paste("\nAbout to try loading package", pkg) )
    
    tryCatch({
      # eval below needed because library() will otherwise be confused
      # https://www.mitchelloharawild.com/blog/loading-r-packages-in-a-loop/
      eval( bquote( library( .(pkg) ) ) )
    }, error = function(err) {
      cat( paste("\n*** COULD NOT LIBRARYIZE PACKAGE:", pkg) )
      any.failed <<- TRUE
    })
    
  }
  if ( any.failed == TRUE ) stop("Some packages couldn't be installed. See outfile for details of which ones.")
  
  # helper code
  path = "/home/groups/manishad/JTE"
  setwd(path)
  source("helper_JTE.R")
  source("stefan_phackR_fns.R")
  
  # ~~ Cluster Run ----------------------------------------
  
  if ( interactive.cluster.run == FALSE ) {
    # get scen parameters (made by genSbatch.R)
    setwd(path)
    scen.params = read.csv( "scen_params.csv" )
    p <<- scen.params[ scen.params$scen == scen, ]
    
    cat("\n\nHEAD OF ENTIRE SCEN.PARAMS:\n")
    print(p)
  }
  
  # ~~ Interactive Cluster Run ----------------------------------------
  # alternatively, generate a simple scen.params in order to run doParallel manually in
  # Sherlock as a test
  if ( interactive.cluster.run == TRUE ) {
    path = "/home/groups/manishad/JTE"
    setwd(path)
    source("helper_JTE.R")
    
    scen.params = tidyr::expand_grid(
      # full list (save):
      #rep.methods = "REML ; ML ; DL ; PMM ; EB ; robu ; jeffreys",
      rep.methods = "jeffreys",
      
      # *If you reorder the args, need to adjust wrangle_agg_local
      ### args shared between sim environments
      k.pub = c(10),  # intentionally out of order so that jobs with boundary choices with complete first
      hack = c("affirm"),
      prob.hacked = c(0),
      # important: if sim.env = stefan, these t2 args are ONLY used for setting start values
      #   and for checking bias of Shat, so set them to have the correct t2a
      #   not clear what t2w should be given the way stefan implements hacking
      #t2a = c(0.05^2, 0.1^2, 0.2^2, 0.5^2, 1),
      t2a = 0.1,
      t2w = c(0),
      # same with Mu
      Mu = c(0, 0.5),
      true.dist = c("norm"),
      
      Nmax = 1,
      m = 50,
      true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"), # original setting  
      rho = c(0),
      
      # Stan control args
      stan.maxtreedepth = 25,
      stan.adapt_delta = 0.995,
      
      get.CIs = TRUE,
      run.optimx = FALSE )
    
    
    scen.params$scen = 1:nrow(scen.params)
    
    scen = 1
  }  # end "if ( interactive.cluster.run == TRUE )"
  
  
  
  # locally, with total k = 100, Nmax = 10, and sim.reps = 250, took 93 min total
  # for that I did sim.reps = 100 per doParallel
  
  # simulation reps to run within this job
  # **this need to match n.reps.in.doParallel in the genSbatch script
  # ***** Set cluster sim reps  -------------------------------------------------
  if ( interactive.cluster.run == FALSE ) sim.reps = 500  # when running all methods except robma
  #if ( interactive.cluster.run == FALSE ) sim.reps = 10  # when running robma only
  
  #if ( interactive.cluster.run == TRUE ) sim.reps = 50 
  
  # set the number of cores
  registerDoParallel(cores=16)
  
}



# FOR LOCAL USE  ------------------------------
if ( run.local == TRUE ) {
  #rm(list=ls())
  
  lapply( toLoad,
          require,
          character.only = TRUE)
  
  
  # helper fns
  code.dir = here()
  setwd(code.dir)
  source("helper_JTE.R")
  
  
  # ~~ ****** Set Local Sim Params -----------------------------
  
  # this is one where Shat behavior was horrible for Jeffreys, but reasonable for other methods
  scen.params = data.frame(
    k.pub = 4,
    t2a = 0.025,
    Mu = 0,
    true.dist = "norm",
    p0 = 0.05,
    Ytype = "cont-SMD",
    N.expr = "40",
    stan.adapt_delta = 0.995,
    stan.maxtreedepth = 25,
    rep.methods = c("MLE-profile ; bayesmeta ; jeffreys")
  )
  
  ### One scen - 105 ###
  # # this is one where Shat behavior was horrible for Jeffreys, but reasonable for other methods
  # scen.params = data.frame(
  #   k.pub = 100,
  #   t2a = 0.0001,
  #   Mu = 0,
  #   true.dist = "norm",
  #   p0 = 0.05,
  #   Ytype = "bin-OR",
  #   N.expr = "40",
  #   stan.adapt_delta = 0.995,
  #   stan.maxtreedepth = 25,
  #   rep.methods = "REML ; DL ; DL2 ; jeffreys"
  # )
  
  # # same, but change dist of SEs
  # # definitely improves jeffreys compared to above, but lower limit is still 0.01, so bad coverage (28%)
  # scen.params = data.frame(
  #   k.pub = 100,
  #   t2a = 0.0001,
  #   Mu = 0,
  #   true.dist = "norm",
  #   p0 = 0.05,
  #   Ytype = "bin-OR",
  #   N.expr = "round( runif(n=1, min=40, max = 4000) )",
  #   stan.adapt_delta = 0.995,
  #  stan.maxtreedepth = 25,
  #   rep.methods = "REML ; DL ; DL2 ; jeffreys"
  # )
  
  # # same, but slightly larger t2a
  # scen.params = data.frame(
  #   k.pub = 100,
  #   t2a = 0.001,
  #   Mu = 0,
  #   true.dist = "norm",
  #   p0 = 0.05,
  #   Ytype = "bin-OR",
  #   N.expr = "round( runif(n=1, min=40, max = 4000) )",
  #   stan.adapt_delta = 0.995,
  #   stan.maxtreedepth = 25,
  #   rep.methods = "REML ; DL ; DL2 ; jeffreys"
  # )
  
  # ### One scen ###
  # scen.params = tidyr::expand_grid(
  #   # full list (save):
  #   #rep.methods = "REML ; ML ; DL ; PMM ; EB ; robu ; jeffreys",
  #   rep.methods = "REML ; DL ; DL2 ; jeffreys",
  #   
  #   # *If you reorder the args, need to adjust wrangle_agg_local
  #   ### args shared between sim environments
  #   k.pub = c(50),  # intentionally out of order so that jobs with boundary choices with complete first
  #  
  #   #t2a = c(0.05^2, 0.1^2, 0.2^2, 0.5^2, 1),
  #   t2a = 0.1,
  # 
  #   # same with Mu
  #   Mu = c(2.3),
  #   true.dist = c("norm"),
  #   
  #   muN = 50,
  #   minN = 20,
  #   Ytype = "bin-OR",
  #   p0 = 0.1,
  #  
  #   # Stan control args
  #   stan.maxtreedepth = 25,
  #   stan.adapt_delta = 0.995,
  #   
  #   get.CIs = TRUE,
  #   run.optimx = FALSE )
  
  # # ### Full set ###
  # ### 2024-01-15 ###
  # scen.params = tidyr::expand_grid(
  #   # full list (save):
  #   rep.methods = "REML ; DL ; DL2 ; PM ; robu ; jeffreys",
  #   
  #   # *If you reorder the args, need to adjust wrangle_agg_local
  #   ### args shared between sim environments
  #   k.pub = c(10,
  #             2, 3, 5, 20, 100),  # intentionally out of order so that jobs with most interesting choices with complete first
  #   
  #   t2a = c(0.1^2, 0.05^2, 0.2^2, 0.5^2),
  #   
  #   # same with Mu
  #   Mu = c(0, 0.5, 1.1, 2.3), # same as Langan's log-ORs
  #   true.dist = c("norm", "expo"),
  #   p0 = c(NA, 0.05, 0.1, 0.5),  
  #   
  #   Ytype = c("cont-SMD", "bin-OR"),
  #   
  #   N.expr = c( "40",
  #               "round( runif(n=1, min=40, max = 400) )",
  #               "400",
  #               "round( runif(n=1, min=2000, max = 4000) )" ),
  #   
  #   # Stan control args
  #   stan.maxtreedepth = 25,
  #   stan.adapt_delta = 0.995,
  #   
  #   get.CIs = TRUE,
  #   run.optimx = FALSE )
  # 
  # table(scen.params$p0, useNA = "ifany")
  # 
  # 
  # #### Remove unwanted combinations
  # 
  # # ... of Mu and Ytype
  # remove = (scen.params$Mu != 0.5) & (scen.params$Ytype == "cont-SMD")
  # scen.params = scen.params[!remove,]
  # # sanity check:
  # table(scen.params$Mu, scen.params$Ytype)
  # 
  # # ... of Ytype and p0
  # remove = rep(FALSE, nrow(scen.params))
  # remove[ !is.na(scen.params$p0) & (scen.params$Ytype == "cont-SMD") ] = TRUE
  # remove[ is.na(scen.params$p0) & (scen.params$Ytype == "bin-OR") ] = TRUE
  # scen.params = scen.params[!remove,]
  # # sanity check:
  # table(scen.params$p0, scen.params$Ytype, useNA = "ifany")
  # #### end of full set
  # 
  
  # add scen numbers
  start.at = 1
  scen.params = scen.params %>% add_column( scen = start.at : ( nrow(scen.params) + (start.at - 1) ),
                                            .before = 1 )
  
  
  ( n.scen = nrow(scen.params) )
  # look at it
  head( as.data.frame(scen.params) )
  ### end of full sim params
  
  scen.params$scen = 1:nrow(scen.params)
  
  # ~ ****** Set local sim.reps  -------------------------------------------------
  sim.reps = 1 
  
  # set the number of local cores
  registerDoParallel(cores=8)
  
  scen = 1
  # data.frame(scen.params %>% filter(scen.name == scen))
  
  # just to avoid errors in doParallel script below
  jobname = "job_1"
  i = 1
}




# COMPILE STAN MODEL ONCE AT BEGINNING------------------------------

if ( run.local == TRUE ) setwd(code.dir)

if ( run.local == FALSE ) setwd(path)


source("init_stan_model_JTE.R")


# RUN SIMULATION ------------------------------

if ( exists("rs") ) rm(rs)

# ~ ****** Beginning of ForEach Loop -----------------------------

# system.time is in seconds
doParallel.seconds = system.time({
  rs = foreach( i = 1:sim.reps, .combine = bind_rows ) %dopar% {
    #for debugging (out file will contain all printed things):
    #for ( i in 1:10 ) {
    
    cat("\n\n~~~~~~~~~~~~~~~~ BEGIN SIM REP", i, "~~~~~~~~~~~~~~~~")
    
    # results for just this simulation rep
    if ( exists("rep.res") ) suppressWarnings( rm(rep.res) )
    
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    cat("\n\n scen variable:\n")
    print(scen)
    
    cat("\n\n scen.params again:\n")
    print(scen.params)
    
    p = scen.params[ scen.params$scen == scen, names(scen.params) != "scen"]
    
    p$V = p$t2a
    p$S = sqrt(p$V)
    
    if ( i == 1 ) cat("\n\nDIM AND HEAD OF P (SINGLE ROW OF SCEN.PARAMS):\n")
    if ( i == 1 ) print(dim(p)); print(p)
    
    # parse methods string
    all.methods = unlist( strsplit( x = p$rep.methods,
                                    split = " ; " ) )
    
    # ~ Simulate Dataset ------------------------------
    # includes unpublished studies
    
    
    d = sim_meta( Mu = p$Mu,
                  t2a = p$t2a,
                  true.dist = p$true.dist,
                  
                  N.expr = p$N.expr,
                  Ytype = p$Ytype,
                  p0 = p$p0,
                  
                  k.pub = p$k.pub)
    
    
    
    d$Zi = d$yi / sqrt(d$vi)
    d$sei = sqrt(d$vi)
    
    if ( i == 1 ) cat("\n\nHEAD OF D:\n")
    if ( i == 1 ) print(head(d))
    
    
    # initialize rep.res st run_method_safe and other standalone estimation fns
    #  will correctly recognize it as having 0 rows
    rep.res = data.frame()
    
    # # ~ Start Values ------------------------------
    Mhat.start = p$Mu
    Shat.start = p$S
    
    # ~ Existing Methods ------------------------------
    
    # ~~ Metafor heterogeneity estimators ------------------------------
    
    # pg 282:
    # https://cran.r-project.org/web/packages/metafor/metafor.pdf
    metafor.methods = all.methods[ all.methods %in% c("REML", "ML", "DL", "EB", "PM", "PMM", "HS", "SJ") ]
    
    
    if ( length(metafor.methods) > 0 ) {
      
      for ( .method in metafor.methods ) {
        rep.res = run_method_safe(method.label = c(.method),
                                  method.fn = function() {
                                    mod = rma( yi = d$yi,
                                               vi = d$vi,
                                               method = .method,
                                               knha = TRUE )
                                    
                                    report_meta(mod, .mod.type = "rma")
                                  },
                                  .rep.res = rep.res )
      }
      
    }
    
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    # ~~ DL2 (Two-Step DL) -------------------------------------------------
    
    # requires 2 calls to metafor, so not handled in the loop above
    if ( "DL2" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("DL2"),
                                method.fn = function() {
                                  
                                  # first step
                                  m0 = rma( yi = d$yi,
                                            vi = d$vi,
                                            method = "DL",
                                            knha = TRUE )
                                  
                                  # second step, per Wolfgang
                                  mod = rma( yi = d$yi,
                                             vi = d$vi,
                                             method = "GENQ",
                                             weights = 1 / (d$vi + m0$tau2),
                                             knha = TRUE )
                                  
                                  report_meta(mod, .mod.type = "rma")
                                },
                                .rep.res = rep.res )
    }
    
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    # ~~ Exact method (package rma.exact) -------------------------------------------------
    
    if ( "exact" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("exact"),
                                method.fn = function() {
                                  
                                  ci = rma.exact.fast( yi = d$yi,
                                                       vi = d$vi )
                                  
                                  # this method doesn't do point estimation of inference for tau
                                  return( list( stats = data.frame( 
                                    MLo = ci[1],
                                    MHi = ci[2]) ) )
                                  
                                  
                                },
                                .rep.res = rep.res )
      
    }
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    # ~~ bayesmeta: Jeffreys prior on *only* tau (package bayesmeta; also Bodnar) -------------------------------------------------
    
    if ( "bayesmeta" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("bayesmeta"),
                                method.fn = function() {
                                  
                                  m = bayesmeta(y = d$yi,
                                                sigma = d$sei,
                                                tau.prior = "Jeffreys")
                                  
                                  # marginal (not joint) intervals
                                  tau_ci = as.numeric( m$post.interval(tau.level=0.95) ) 
                                  mu_ci = as.numeric( m$post.interval(mu.level=0.95) )
                                  
                                  # this method doesn't do point estimation of inference for tau
                                  return( list( stats = data.frame( 
                                    Mhat = m$MAP["joint", "mu"],
                                    Shat = m$MAP["joint", "tau"],
                                    MLo = mu_ci[1],
                                    MHi = mu_ci[2],
                                    SLo = tau_ci[1],
                                    SHi = tau_ci[2] ) ) )
                                  
                                  
                                },
                                .rep.res = rep.res )
      
    }
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    # ~~ Own stan model: Jeffreys prior on *only* tau ------------------------------
    
    if ( "jeffreys-tau" %in% all.methods ) {
      # this one has multiple labels in method arg because a single call to estimate_jeffreys_mcmc
      #  returns 2 lines of output, one for posterior mean and one for posterior median
      # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
      rep.res = run_method_safe(method.label = c("jeffreys-tau-pmean",
                                                 "jeffreys-tau-pmed",
                                                 "jeffreys-tau-max-lp-iterate"),
                                method.fn = function() estimate_jeffreys_tau_only(.yi = d$yi,
                                                                                  .sei = d$sei,
                                                                                  
                                                                                  .Mu.start = Mhat.start,
                                                                                  # can't handle start value of 0:
                                                                                  .Tt.start = max(0.01, Shat.start),
                                                                                  .stan.adapt_delta = p$stan.adapt_delta,
                                                                                  .stan.maxtreedepth = p$stan.maxtreedepth), .rep.res = rep.res )
      
      
      # start values for finding posterior mode analytically
      maxlp = rep.res[ rep.res$method == "jeffreys-tau-max-lp-iterate", ]
      Mhat.MaxLP = maxlp$Mhat
      Shat.MaxLP = maxlp$Shat
      
      
      # find posterior mode analytically
      rep.res = run_method_safe(method.label = c("jeffreys-tau-pmode"),
                                method.fn = function() {
                                  
                                  # as in the future pkg
                                  mle_fit <- mle_params(Mhat.MaxLP, Shat.MaxLP, d$yi, d$sei)
                                  modes <- c(mle_fit@coef[["mu"]], mle_fit@coef[["tau"]])
                                  optim_converged <- mle_fit@details$convergence == 0
                                  
                                  
                                  return( list( stats = data.frame( 
                                    
                                    Mhat = modes[1],
                                    Shat = modes[2],
                                    
                                    # all inference is again from MCMC
                                    MhatSE = maxlp$MhatSE,
                                    ShatSE = maxlp$ShatSE,
                                    MLo = maxlp$MLo,
                                    MHi =maxlp$MHi,
                                    SLo = maxlp$SLo,
                                    SHi = maxlp$SHi,
                                    
                                    stan.warned = maxlp$stan.warned,
                                    stan.warning = maxlp$stan.warning,
                                    MhatRhat = maxlp$MhatRhat,
                                    ShatRhat = maxlp$ShatRhat,
                                    
                                    OptimConverged = optim_converged) ) )
                                  
                                }, .rep.res = rep.res )
      
      
      
      cat("\n doParallel flag: Done jeffreys-tau if applicable")
    }
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    
    # ~~ MLE-profile ------------------------------
    
    if ( "MLE-profile" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("MLE-profile"),
                                method.fn = function() {
                                  
                                  nll_fun <- function(mu, tau) get_nll(mu, tau, d$yi, d$sei)
                                  my_mle = stats4::mle(minuslogl = nll_fun,
                                                       start = list(mu = Mhat.start, tau = Shat.start),
                                                       method = "L-BFGS-B")
                                  
                                  # this fn will complain if optimizer in my_mle hasn't found the true max
                                  #  which is a good thing
                                  cis = stats4::confint(my_mle)
                                  
                                  return( list( stats = data.frame( 
                                    
                                    Mhat = as.numeric( attr(my_mle, "coef")["mu"] ),
                                    Shat = as.numeric( attr(my_mle, "coef")["tau"] ),
                                    
                                    MhatSE = NA,
                                    ShatSE = NA,
                                    
                                    MLo = cis["mu", 1],
                                    MHi = cis["mu", 2],
                                    
                                    SLo = max( 0, cis["tau", 1] ),
                                    SHi = cis["tau", 2] ) ) )
                                  
                                  
                                },
                                .rep.res = rep.res )
      
    }
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    # ~~ Robust Variance Estimation ------------------------------
    
    # when no clustering or userweights, coincides exactly with DL
    #  so only useful as a sanity check
    if ( "robu" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("robu"),
                                method.fn = function() {
                                  mod = robu( yi ~ 1,
                                              data = d,
                                              studynum = 1:nrow(d),
                                              var.eff.size = vi,
                                              small = TRUE)
                                  
                                  report_meta(mod, .mod.type = "robu")
                                },
                                .rep.res = rep.res )
    }
    
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    # ~ New Methods ------------------------------
    
    # ~~ ***** Jeffreys (MCMC) ------------------------------
    
    if ( "jeffreys" %in% all.methods ) {
      # this one has multiple labels in method arg because a single call to estimate_jeffreys_mcmc
      #  returns 2 lines of output, one for posterior mean and one for posterior median
      # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
      rep.res = run_method_safe(method.label = c("jeffreys-pmean",
                                                 "jeffreys-pmed",
                                                 "jeffreys-max-lp-iterate"),
                                method.fn = function() estimate_jeffreys(.yi = d$yi,
                                                                         .sei = d$sei,
                                                                         
                                                                         .Mu.start = Mhat.start,
                                                                         # can't handle start value of 0:
                                                                         .Tt.start = max(0.01, Shat.start),
                                                                         .stan.adapt_delta = p$stan.adapt_delta,
                                                                         .stan.maxtreedepth = p$stan.maxtreedepth), .rep.res = rep.res )
      
      
      # start values for finding posterior mode analytically
      maxlp = rep.res[ rep.res$method == "jeffreys-max-lp-iterate", ]
      Mhat.MaxLP = maxlp$Mhat
      Shat.MaxLP = maxlp$Shat
      
      
      # find posterior mode analytically
      rep.res = run_method_safe(method.label = c("jeffreys-pmode"),
                                method.fn = function() {
                                  
                                  # as in the future pkg
                                  mle_fit <- mle_params(Mhat.MaxLP, Shat.MaxLP, d$yi, d$sei)
                                  modes <- c(mle_fit@coef[["mu"]], mle_fit@coef[["tau"]])
                                  optim_converged <- mle_fit@details$convergence == 0
                                  
                                  
                                  return( list( stats = data.frame( 
                                    
                                    Mhat = modes[1],
                                    Shat = modes[2],
                                    
                                    # all inference is again from MCMC
                                    MhatSE = maxlp$MhatSE,
                                    ShatSE = maxlp$ShatSE,
                                    MLo = maxlp$MLo,
                                    MHi =maxlp$MHi,
                                    SLo = maxlp$SLo,
                                    SHi = maxlp$SHi,
                                    
                                    stan.warned = maxlp$stan.warned,
                                    stan.warning = maxlp$stan.warning,
                                    MhatRhat = maxlp$MhatRhat,
                                    ShatRhat = maxlp$ShatRhat,
                                    
                                    OptimConverged = optim_converged) ) )
                                  
                                }, .rep.res = rep.res )
      
      
      
      cat("\n doParallel flag: Done jeffreys if applicable")
    }
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    # ~ Add Scen Params and Sanity Checks  -------------------------------------------------
    
    # add in scenario parameters
    # do NOT use rbind here; bind_cols accommodates possibility that some methods' rep.res
    #  have more columns than others
    rep.res = p %>% bind_cols( rep.res )
    
    # add more info
    rep.res = rep.res %>% add_column( rep.name = i, .before = 1 )
    rep.res = rep.res %>% add_column( scen.name = scen, .before = 1 )
    rep.res = rep.res %>% add_column( job.name = jobname, .before = 1 )
    
    cat("\ndoParallel flag: Before adding sanity checks to rep.res")
    rep.res$sancheck_mean_yi = mean(d$mean_yi)
    rep.res$sancheck_mean_pY0 = mean(d$pY0)
    rep.res$sancheck_mean_pY1 = mean(d$pY1)
    rep.res$sancheck_mean_pY = mean(d$pY)
    
    rep.res$sancheck_mean_nY0 = mean(d$nY0)
    rep.res$sancheck_mean_nY0_theory = mean(d$nY0_theory)
    rep.res$sancheck_mean_nY1 = mean(d$nY1)
    rep.res$sancheck_mean_nY1_theory = mean(d$nY1_theory)
    
    rep.res$sancheck_mean_N = mean(d$N)
    
    rep.res
    
  }  ### end foreach loop
  
} )[3]  # end system.time


# quick look
#rs %>% dplyr::select(method, Shat, SLo, SHi, Mhat, MLo, MHi)

names(rs)
table(rs$method)



# ~~ End of ForEach Loop ----------------
# estimated time for 1 simulation rep
# use NAs for additional methods so that the SUM of the rep times will be the
#  total computational time
nMethods = length( unique(rs$method) )

print(nMethods)
print(unique(rs$method))
table(rs$method)

print(nrow(rs))

rs$doParallel.seconds = doParallel.seconds


rs$rep.seconds = doParallel.seconds/sim.reps
rs$rep.seconds[ rs$method != unique(rs$method)[1] ] = NA

expect_equal( as.numeric( sum(rs$rep.seconds, na.rm = TRUE) ),
              as.numeric(doParallel.seconds) )


# ~ QUICK RESULTS SUMMARY ---------------------------------------------------

if ( run.local == TRUE ) {
  
  Mu = p$Mu
  tau = sqrt(p$t2a)
  
  # quick look locally
  # rs %>% mutate_if(is.numeric, function(x) round(x,2) )
  
  agg = rs %>% group_by(method) %>%
    summarise( PropMhatNA = mean(is.na(Mhat)),
               PropCI.NA = mean(is.na(MLo)),
               
               Shat = meanNA(Shat),
               ShatMSE = meanNA( (Shat - tau)^2 ),
               ShatBias = meanNA(Shat - tau),
               ShatSE = meanNA(ShatSE),
               #ShatMn = meanNA(Shat),
               
               ShatCover = meanNA(SLo <= tau & SHi >= tau),
               ShatWidth = meanNA( SHi - SLo ),
               SLo = meanNA(SLo),
               SHi = meanNA(SHi),
               
               Mhat = meanNA(Mhat),
               MhatMSE = meanNA( (Mhat - Mu)^2 ),
               MhatBias = meanNA(Mhat - Mu),
               MhatEmpSE = sd( Mhat, na.rm = TRUE ),
               #ShatMn = meanNA(Shat),
               
               MhatCover = meanNA(MLo <= Mu & MHi >= Mu),
               MhatWidth = meanNA( MHi - MLo ),
               MLo = meanNA(MLo),
               MHi = meanNA(MHi),
               
               sancheck_mean_yi = meanNA(sancheck_mean_yi),
               sancheck_mean_pY0 = meanNA(sancheck_mean_pY0),
               sancheck_mean_pY1 = meanNA(sancheck_mean_pY1),
               
               sancheck_mean_nY0 = meanNA(sancheck_mean_nY0),
               sancheck_mean_nY0_theory = meanNA(sancheck_mean_nY0_theory),
               
               sancheck_mean_nY1 = meanNA(sancheck_mean_nY1),
               sancheck_mean_nY1_theory = meanNA(sancheck_mean_nY1_theory),
               mean_N = mean(sancheck_mean_N) )
  
  # round
  agg = as.data.frame( agg %>% mutate_if( is.numeric,
                                          function(x) round(x,2) ) )
  
  agg
  
  
  # scenario diagnostics for scenario
  keepers = namesWith("sancheck_", rs)
  agg.checks = rs %>% summarise_at( keepers,
                                    function(x) round( mean(x), 2) )
  
  t(agg.checks)
  
}




# ~ WRITE LONG RESULTS ------------------------------
if ( run.local == FALSE ) {
  setwd("/home/groups/manishad/JTE/long_results")
  fwrite( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
}