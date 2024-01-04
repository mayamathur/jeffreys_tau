
# IMPORTANT NOTES -----------------------------

# Important things to remember: 
#
# - The returned Vhat is an estimate of T2 + t2w, *not* T2 itself

# Debugging help:
# 
# - The jobs may fail before fitting modAll with no apparent errors if 
#   k is too large for rma.mv. In that case, try setting p$k < 500 for modAll
#  and modPub small to prevent those models from being fit. 


# for interactive Sherlock:
# path = "/home/groups/manishad/SAPH"
# setwd(path)
# source("doParallel_SAPH.R")


# because Sherlock 2.0 restores previous workspace
rm( list = ls() )


# are we running locally?
run.local = FALSE

# should we set scen params interactively on cluster?
interactive.cluster.run = FALSE

# ~~ Packages -----------------------------------------------
toLoad = c("crayon",
           "dplyr",
           "foreach",
           "doParallel",
           "boot",
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
           "truncdist",
           "tibble",
           "tmvtnorm",
           "testthat",
           "truncreg",
           "truncnorm",
           "rstan", # note: to reinstall this one, need to use high-mem session
           "optimx",
           "weightr",
           "phacking",
           "RoBMA")  # note: to reinstall this one, need ml load jags

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
  path = "/home/groups/manishad/SAPH"
  setwd(path)
  source("helper_SAPH.R")
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
    path = "/home/groups/manishad/SAPH"
    setwd(path)
    source("helper_SAPH.R")
    source("stefan_phackR_fns.R")
    scen.params = tidyr::expand_grid(
      rep.methods = "naive ; jeffreys-mcmc ; jeffreys-sd",
      
      # args from sim_meta_2
      Nmax = 30,
      Mu = c(0.5),
      t2a = 0,
      t2w = 0,
      m = 50,
      
      
      hack = c("favor-best-affirm-wch"),
      rho = c(0),
      k.pub.nonaffirm = 20,
      prob.hacked = c(0.8),
      
      true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
      
      # Stan control args
      stan.maxtreedepth = 20,
      stan.adapt_delta = 0.98,
      
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
  if ( interactive.cluster.run == FALSE ) sim.reps = 25  # when running all methods except robma
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
  code.dir = here("Sherlock code")
  setwd(code.dir)
  source("helper_SAPH.R")
  source("stefan_phackR_fns.R")
  
  
  # ~~ ****** Set Local Sim Params -----------------------------
  
  # # Mathur environment
  # scen.params = tidyr::expand_grid(
  #   # full list (save):
  #   #rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; jeffreys-mcmc ; jeffreys-sd ; prereg-naive",
  #   rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; pet-peese ; robma ; jeffreys-mcmc ; prereg-naive",
  #   #rep.methods = "naive",
  #   #rep.methods = "naive ; gold-std ; jeffreys-mcmc ; 2psm",
  #   #rep.methods = c("rtma-pkg ; jeffreys-mcmc"),
  #   #rep.methods = "naive ; rtma-pkg",
  #   
  #   sim.env = "mathur",
  #   
  #   # *If you reorder the args, need to adjust wrangle_agg_local
  #   ### args shared between sim environments
  #   k.pub.nonaffirm = c(10),  # intentionally out of order so that jobs with boundary choices with complete first 
  #   hack = c("affirm"),
  #   prob.hacked = c(0.8),
  #   # important: if sim.env = stefan, these t2 args are ONLY used for setting start values
  #   #   and for checking bias of Shat, so set them to have the correct t2a
  #   #   not clear what t2w should be given the way stefan implements hacking 
  #   t2a = c(0),
  #   t2w = c(0.2^2),
  #   # same with Mu
  #   Mu = c(0.5),
  #   
  #   ### only needed if sim.env = "mathur": args from sim_meta_2
  #   Nmax = 30,
  #   m = 50,
  #   true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
  #   rho = c(0),
  #   ### end of stuff for sim.env = "mathur"
  #   
  #   # ### only needed if sim.env = "stefan": args from sim_meta_2
  #   # strategy.stefan = c("firstsig", "smallest"),  # "firstsig" or "smallest"
  #   # alternative.stefan = c("greater", "two.sided"),  # "two.sided" or "greater"
  #   # stringent.hack = TRUE,  # mathur sims always effectively use stringent.hack = TRUE
  #   # ### end of stuff for sim.env = "stefan"
  #   
  #   # Stan control args
  #   stan.maxtreedepth = 25,
  #   stan.adapt_delta = 0.995,
  #   
  #   get.CIs = TRUE,
  #   run.optimx = FALSE )
  
  
  
  # Stefan environment - exactly as in genSbatch
  scen.params = tidyr::expand_grid(
    # without robma:
    rep.methods = "naive ; gold-std ; maon",
    #rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; pet-peese ; rtma-pkg ; prereg-naive",
    
    sim.env = "stefan",
    
    ### args shared between sim environments
    #k.pub.nonaffirm = c(10, 15, 20, 30, 50, 70, 100),
    k.pub.nonaffirm = c(100),  # intentionally out of order so that jobs with boundary choices with complete first
    hack = c("DV"),
    prob.hacked = c(0.8),
    # important: if sim.env = stefan, these t2 args are ONLY used for setting start values
    #   and for checking bias of Shat, so set them to have the correct t2a
    #   not clear what t2w should be given the way stefan implements hacking
    t2a = c(1),
    t2w = c(0),
    # same with Mu
    Mu = c(0),
    
    # ### only needed if sim.env = "mathur": args from sim_meta_2
    # Nmax = 30,
    # m = 50,
    #
    # true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
    # rho = c(0),
    # ### end of stuff for sim.env = "mathur"
    
    ### only needed if sim.env = "stefan": args from sim_meta_2
    strategy.stefan = c("firstsig"),  # "firstsig" or "smallest"
    alternative.stefan = c("greater"),  # "two.sided" or "greater"
    stringent.hack = TRUE,  # mathur sims always effectively use stringent.hack = TRUE
    ### end of stuff for sim.env = "stefan"
    
    # Stan control args
    stan.maxtreedepth = 25,
    stan.adapt_delta = 0.995,
    
    get.CIs = TRUE,
    run.optimx = FALSE )
  
  
  
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

# READ IN LODDER SEs ONCE AT BEGINNING ------------------------------

# only needed if using "draw_lodder_se()" as one of the true.sei.expr

if ( "draw_lodder_se()" %in% scen.params$true.sei.expr ) {
  
  setwd("/home/groups/manishad/SAPH/applied_examples/data")
  d.lodder = fread("lodder_prepped.csv")
  lodder.ses = sqrt(d.lodder$vi)
  
  cat("\n\ndoParallel: just read in Lodder SEs:")
  summary(lodder.ses)
  
}



# COMPILE STAN MODEL ONCE AT BEGINNING------------------------------

if ( run.local == TRUE ) setwd(code.dir)

if ( run.local == FALSE ) setwd(path)


source("init_stan_model_SAPH.R")


# RUN SIMULATION ------------------------------

if ( exists("rs") ) rm(rs)

# ~ Beginning of ForEach Loop -----------------------------

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
    
    # calculate TOTAL heterogeneity
    p$V = p$t2a + p$t2w
    p$S = sqrt(p$V)
    
    if ( i == 1 ) cat("\n\nDIM AND HEAD OF P (SINGLE ROW OF SCEN.PARAMS):\n")
    if ( i == 1 ) print(dim(p)); print(p)
    
    # parse methods string
    all.methods = unlist( strsplit( x = p$rep.methods,
                                    split = " ; " ) )
    
    # ~ Simulate Dataset ------------------------------
    # includes unpublished studies
    
    if ( p$sim.env == "mathur" ) {
      d = sim_meta_2( Nmax = p$Nmax,
                      Mu = p$Mu,
                      t2a = p$t2a,
                      m = p$m,
                      t2w = p$t2w,
                      true.sei.expr = p$true.sei.expr,
                      hack = p$hack,
                      rho = p$rho,
                      
                      k.pub.nonaffirm = p$k.pub.nonaffirm,
                      prob.hacked = p$prob.hacked,
                      return.only.published = FALSE)
    }
    
    
    if ( p$sim.env == "stefan" ) {
      d = sim_meta_2_stefan( hack.type = p$hack,
                             stringent.hack = p$stringent.hack,
                             strategy.stefan = p$strategy.stefan,
                             alternative.stefan = p$alternative.stefan,
                             
                             k.pub.nonaffirm = p$k.pub.nonaffirm,
                             prob.hacked = p$prob.hacked,
                             return.only.published = FALSE)
      
      # recode for comparability with sim_meta_2
      d$hack = NA
      d$hack[ d$is.hacked == 1 ] = p$hack
      d$hack[ d$is.hacked == 0 ] = "no"
      
      
      d %>% group_by(Di, is.hacked, affirm) %>%
        summarise( n(),
                   mean(yi),
                   mean(vi) )
      
      hist(d$yi)
      
      # # sanity checks
      # hist(d$yi[ d$Di == 1] )
      # hist(d$yi[ d$Di == 1 & d$is.hacked == TRUE ] )
      # hist(d$yi[ d$Di == 1 & d$is.hacked == TRUE ]/ d$sei[ d$Di == 1 & d$is.hacked == TRUE ] ) # Z-scores
      # 
      # hist(d$yi[ d$Di == 1 & d$is.hacked == FALSE ] )
      # hist(d$yi[ d$Di == 1 & d$is.hacked == FALSE ]/ d$sei[ d$Di == 1 & d$is.hacked == FALSE ] ) # Z-scores
      
    }
    
    
    d$Zi = d$yi / sqrt(d$vi)
    
    # ~ Dataset Subsets for Various Methods ------------------------------
    # dataset of only published results
    dp = d %>% filter(Di == 1)
    dp$sei = sqrt(dp$vi)
    
    # keep first draws only
    # for Stefan, this is the same as d
    d.first = d[ !duplicated(d$study), ]
    
    # published nonaffirmatives only
    dpn = dp[ dp$affirm == FALSE, ]
    
    # published affirmatives only
    dpa = dp[ dp$affirm == TRUE, ]
    
    # discard affirmatives from hacked studies
    dp.csm = dp %>% filter( affirm == FALSE | hack == "no" )
    
    # this is like analyzing only preregistered studies
    dp.unhacked = dp %>% filter(hack == "no")
    
    if ( i == 1 ) cat("\n\nHEAD OF DP:\n")
    if ( i == 1 ) print(head(dp))
    
    # # published hacked ones only
    # dph = d %>% filter(hack == p$hack & Di == 1)
    # 
    # # published unhacked ones only
    # # so only one per study set
    # # same as second row of above table
    # dpu = d %>% filter(hack == "no" & Di == 1)
    
    # initialize rep.res st run_method_safe and other standalone estimation fns
    #  will correctly recognize it as having 0 rows
    rep.res = data.frame()
    
    # ~ Start Values ------------------------------
    Mhat.start = p$Mu
    Shat.start = p$S
    
    # in case we're not doing jeffreys-mcmc or it fails
    Mhat.MaxLP = NA
    Shat.MaxLP = NA
    
    Mhat.MAP = NA
    Shat.MAP = NA
    
    
    # ~ Existing Methods ------------------------------
    
    # ~~ Naive Meta-Analysis (All PUBLISHED Draws) ------------------------------
    
    if ( "naive" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("naive"),
                                method.fn = function() {
                                  mod = rma( yi = dp$yi,
                                             vi = dp$vi,
                                             method = "REML",
                                             knha = TRUE )
                                  
                                  report_meta(mod, .mod.type = "rma")
                                },
                                .rep.res = rep.res )
    }
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
    # ~~ Gold-Standard Meta-Analysis (ALL FIRST Draws) ------------------------------
    
    if ( "gold-std" %in% all.methods ) {
      
      if ( p$sim.env == "mathur" ) {
        rep.res = run_method_safe(method.label = c("gold-std"),
                                  method.fn = function() {
                                    mod.all = rma( yi = d.first$yi,
                                                   vi = d.first$vi,
                                                   method = "REML",
                                                   knha = TRUE )
                                    
                                    report_meta(mod.all, .mod.type = "rma")
                                  },
                                  .rep.res = rep.res )
      }
      
      # Stefan datasets are wide (each row has yi [hacked] and yio [ideal draw])
      #  rather than long, so handled slightly differently
      if ( p$sim.env == "stefan" ) {
        rep.res = run_method_safe(method.label = c("gold-std"),
                                  method.fn = function() {
                                    mod.all = rma( yi = d.first$yio,
                                                   vi = d.first$seio^2,
                                                   method = "REML",
                                                   knha = TRUE )
                                    
                                    report_meta(mod.all, .mod.type = "rma")
                                  },
                                  .rep.res = rep.res )
      }
    }
    
    
    if (run.local == TRUE) srr(rep.res)
    
    # ~~ MAON (Nonaffirmative Published Draws) ------------------------------
    
    if ( "maon" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("maon"),
                                method.fn = function() {
                                  mod = robu( yi ~ 1, 
                                              data = dpn, 
                                              studynum = 1:nrow(dpn),
                                              var.eff.size = vi,
                                              small = TRUE )
                                  
                                  report_meta(mod, .mod.type = "robu")
                                },
                                .rep.res = rep.res )
      
    }
    
    
    # ~~ 2PSM (All Published Draws) ------------------------------
    
    if ( "2psm" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("2psm"),
                                method.fn = function() {
                                  mod = weightfunct( effect = dp$yi,
                                                     v = dp$vi,
                                                     steps = c(0.025, 1),
                                                     table = TRUE ) 
                                  
                                  H = mod[[2]]$hessian
                                  ses = sqrt( diag( solve(H) ) )
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = mod[[2]]$par[2],
                                                            MLo = mod[[2]]$par[2] - qnorm(.975) * ses[2],
                                                            MHi = mod[[2]]$par[2] + qnorm(.975) * ses[2],
                                                            
                                                            Shat = sqrt( mod[[2]]$par[1] ),
                                                            # truncate lower limit at 0
                                                            SLo = sqrt( max( 0, mod[[2]]$par[1] - qnorm(.975) * ses[1] ) ),
                                                            SHi = sqrt( mod[[2]]$par[1] + qnorm(.975) * ses[1] ) ) ) 
                                },
                                .rep.res = rep.res )
      
    }
    
    
    # ~~ P-Curve (Published Affirmatives) ------------------------------
    
    if ( "pcurve" %in% all.methods ) {
      # since pcurve.opt internally retains only 
      rep.res = run_method_safe(method.label = c("pcurve"),
                                method.fn = function() {
                                  #later, revisit the decision to use df_obs = 1000
                                  #   to effectively treat yi/sei z-scores
                                  
                                  # using all significant studies (that's what pcurve.opt does internally):
                                  Mhat = pcurve.opt( t_obs = dp$Zi,
                                                     df_obs = rep(1000, length(dp$Zi)),
                                                     dmin = -5, #HARD-CODED and arbitrary
                                                     dmax = 5)
                                  
                                  
                                  # # using only published affirmatives:
                                  # Mhat = pcurve.opt( t_obs = dpa$yi/dpa$sei,
                                  #                    df_obs = rep(1000, length(dpa$yi)),
                                  #                    dmin = -5, #HARD-CODED and arbitrary
                                  #                    dmax = 5)
                                  
                                  
                                  return( list( stats = data.frame( Mhat = Mhat) ) )
                                },
                                .rep.res = rep.res )
      
    }
    
    
    
    # ~~ PET-PEESE ------------------------------
    
    # best PET-PEESE reference: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pub-bias.html
    if ( "pet-peese" %in% all.methods ) {
      # c.f. Hilgard's repo: https://github.com/Joe-Hilgard/PETPEESE/blob/master/PETPEESE_functions.R
      rep.res = run_method_safe(method.label = c("pet-peese"),
                                method.fn = function() {
                                  
                                  # first fit PET
                                  mod.PET = lm(yi ~ sei,
                                               weights = 1/vi,
                                               data = dp)
                                  
                                  # check if coefficient for intercept is significant
                                  # needs to be a one-sided test
                                  PET.pval.twotail = summary(mod.PET)$coefficients["(Intercept)", "Pr(>|t|)"]
                                  PET.intercept = summary(mod.PET)$coefficients["(Intercept)", "Estimate"]
                                  
                                  
                                  # if PET is significant, run PEESE instead
                                  # alpha = 0.10 for two-tailed p-value to have one-tailed alpha = 0.05
                                  if ( PET.pval.twotail < 0.10 & PET.intercept > 0 ) {
                                    mod.PEESE = lm(yi ~ vi, # regress on vi instead of sei 
                                                   weights = 1/vi,
                                                   data = dp)
                                    
                                    mod = mod.PEESE
                                  } else {
                                    # in this case, stick with PET
                                    mod = mod.PET
                                  }
                                  
                                  # follow the same return structure as report_meta
                                  summ = summary(mod)
                                  CIs = confint(mod)
                                  list( stats = data.frame( Mhat = as.numeric( mod$coefficients[1] ),
                                                            MLo = as.numeric( CIs[1,1] ),
                                                            MHi = as.numeric( CIs[1,2] ) ) )
                                },
                                .rep.res = rep.res ) 
      
    }
    
    
    
    # ~~ RoBMA ------------------------------
    # this method is very slow! 
    
    if (run.local == TRUE) srr(rep.res)
    
    # https://cran.r-project.org/web/packages/RoBMA/vignettes/CustomEnsembles.html
    if ( "robma" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("robma"),
                                method.fn = function() {
                                  
                                  # default model ensemble
                                  mod = RoBMA(y = dp$yi,
                                              se = dp$sei)
                                  
                                  
                                  summ  = summary(mod)
                                  
                                  # does also return tau, but I'm not including it for now
                                  list( stats = data.frame( Mhat = as.numeric( mod$coefficients[1] ),
                                                            MLo = as.numeric( summ$estimates["mu", "0.025"] ),
                                                            MHi = as.numeric( summ$estimates["mu", "0.975"] ) ) )
                                },
                                .rep.res = rep.res ) 
      
    }
    
    if (run.local == TRUE) srr(rep.res)
    
    # ~ New Methods ------------------------------
    # ~~ ***** MCMC ------------------------------
    
    if ( "jeffreys-mcmc" %in% all.methods ) {
      # # temp for refreshing code
      # path = "/home/groups/manishad/SAPH"
      # setwd(path)
      # source("helper_SAPH.R")
      # source("init_stan_model_SAPH.R")
      # 
      # 
      # #TEMP
      # estimate_jeffreys_mcmc_RTMA(.yi = dpn$yi,
      #                             .sei = sqrt(dpn$vi),
      #                             .tcrit = dpn$tcrit,
      #                             .Mu.start = Mhat.start,
      #                             # can't handle start value of 0:
      #                             .Tt.start = max(0.01, Shat.start),
      #                             .stan.adapt_delta = p$stan.adapt_delta,
      #                             .stan.maxtreedepth = p$stan.maxtreedepth)
      
      # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
      #  returns 2 lines of output, one for posterior mean and one for posterior median
      # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
      rep.res = run_method_safe(method.label = c("jeffreys-mcmc-pmean",
                                                 "jeffreys-mcmc-pmed",
                                                 "jeffreys-mcmc-max-lp-iterate"),
                                method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dpn$yi,
                                                                                   .sei = sqrt(dpn$vi),
                                                                                   .tcrit = dpn$tcrit,
                                                                                   .Mu.start = Mhat.start,
                                                                                   # can't handle start value of 0:
                                                                                   .Tt.start = max(0.01, Shat.start),
                                                                                   .stan.adapt_delta = p$stan.adapt_delta,
                                                                                   .stan.maxtreedepth = p$stan.maxtreedepth), .rep.res = rep.res )
      
      
      Mhat.MaxLP = rep.res$Mhat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
      Shat.MaxLP = rep.res$Shat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
      
      cat("\n doParallel flag: Done jeffreys-mcmc if applicable")
    }
    
    
    # ~~ ***** RTMA using package ------------------------------
    
    # uses posterior mode
    if ( "rtma-pkg" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("rtma-pkg"),
                                method.fn = function() {
                                  
                                  # important: when running stefan sim env, need to subset on affirm status here because they sometimes use one-tailed selection
                                  mod = phacking_meta( yi = dp$yi[dp$affirm == FALSE],
                                                       vi = dp$vi[dp$affirm == FALSE],
                                                       stan_control = list(adapt_delta = p$stan.adapt_delta,
                                                                           max_treedepth = p$stan.maxtreedepth),
                                                       parallelize = FALSE )
                                  
                                  stats = mod$stats
                                  
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = stats$mode[stats$param == "mu"],
                                                            Shat =  stats$mode[stats$param == "tau"],
                                                            
                                                            MhatSE = stats$se[stats$param == "mu"],
                                                            ShatSE = stats$se[stats$param == "tau"],
                                                            
                                                            MLo = stats$ci_lower[ stats$param == "mu"],
                                                            MHi = stats$ci_upper[ stats$param == "mu"],
                                                            
                                                            
                                                            SLo =  stats$ci_lower[stats$param == "tau"],
                                                            SHi = stats$ci_upper[ stats$param == "tau"],
                                                            
                                                            MhatRhat = stats$r_hat[stats$param == "mu"],
                                                            ShatRhat = stats$r_hat[stats$param == "tau"]
                                  ) ) 
                                  
                                },
                                .rep.res = rep.res )
      
      
      cat("\n doParallel flag: Done jeffreys-mcmc if applicable")
    }
    
    
    if ( run.local == TRUE ) srr(rep.res)
    
    # ~~ Change Starting Values -----
    if ( !is.na(Mhat.MaxLP) ) Mhat.start = Mhat.MaxLP 
    if ( !is.na(Shat.MaxLP) ) Shat.start = Shat.MaxLP 
    
    # ~~ MAP (SD param) ------------------------------
    
    #bm
    if ( "jeffreys-sd" %in% all.methods ) {
      # # temp for refreshing code
      # path = "/home/groups/manishad/SAPH"
      # setwd(path)
      # source("helper_SAPH.R")
      # source("init_stan_model_SAPH.R")
      
      rep.res = run_method_safe(method.label = c("jeffreys-sd"),
                                method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                              sei = sqrt(dpn$vi),
                                                                              par2is = "Tt",
                                                                              tcrit = dpn$tcrit, 
                                                                              Mu.start = Mhat.start,
                                                                              par2.start = Shat.start,
                                                                              usePrior = TRUE,
                                                                              get.CIs = p$get.CIs,
                                                                              CI.method = "wald",
                                                                              run.optimx = p$run.optimx
                                ),
                                .rep.res = rep.res )
      
      Mhat.MAP = rep.res$Mhat[ rep.res$method == "jeffreys-sd" ]
      Shat.MAP = rep.res$Shat[ rep.res$method == "jeffreys-sd" ]
    }
    
    
    # ~~ MAP (Var param) ------------------------------
    
    if ( "jeffreys-var" %in% all.methods ) {
      # # temp for refreshing code
      # path = "/home/groups/manishad/SAPH"
      # setwd(path)
      # source("helper_SAPH.R")
      
      rep.res = run_method_safe(method.label = c("jeffreys-var"),
                                method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                              sei = sqrt(dpn$vi),
                                                                              par2is = "T2t",
                                                                              tcrit = dpn$tcrit, 
                                                                              Mu.start = Mhat.start,
                                                                              par2.start = Shat.start^2,
                                                                              usePrior = TRUE,
                                                                              get.CIs = p$get.CIs,
                                                                              CI.method = "wald",
                                                                              
                                                                              run.optimx = p$run.optimx),
                                .rep.res = rep.res )
      
      Mhat.MAP = rep.res$Mhat[ rep.res$method == "jeffreys-var" ]
      Shat.MAP = rep.res$Shat[ rep.res$method == "jeffreys-var" ]
    }
    
    # ~~ Change Starting Values -----
    
    if ( !is.na(Mhat.MAP) ) Mhat.start = Mhat.MAP 
    if ( !is.na(Shat.MAP) ) Shat.start = Shat.MAP 
    
    # ~~ MLE (SD param) ------------------------------
    
    if ( "mle-sd" %in% all.methods ) {
      
      # # temp for refreshing code
      # path = "/home/groups/manishad/SAPH"
      # setwd(path)
      # source("helper_SAPH.R")
      
      rep.res = run_method_safe(method.label = c("mle-sd"),
                                method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                              sei = sqrt(dpn$vi),
                                                                              par2is = "Tt",
                                                                              tcrit = dpn$tcrit, 
                                                                              Mu.start = Mhat.start,
                                                                              par2.start = Shat.start,
                                                                              usePrior = FALSE,
                                                                              get.CIs = p$get.CIs,
                                                                              CI.method = "wald",
                                                                              run.optimx = p$run.optimx),
                                .rep.res = rep.res )
      
      
      
      cat("\n doParallel flag: Done mle-sd if applicable")
    }
    
    # ~~ CSM: MLE with known hacking status (SD param) ------------------------------
    
    
    
    if ( "csm-mle-sd" %in% all.methods &
         nrow(dp.csm) > 0 ) {
      
      # # temp for refreshing code
      # path = "/home/groups/manishad/SAPH"
      # setwd(path)
      # source("helper_SAPH.R")
      
      rep.res = run_method_safe(method.label = c("csm-mle-sd"),
                                method.fn = function() estimate_jeffreys_RTMA(yi = dp.csm$yi,
                                                                              sei = sqrt(dp.csm$vi),
                                                                              par2is = "Tt",
                                                                              tcrit = qnorm(0.975), 
                                                                              Mu.start = Mhat.start,
                                                                              par2.start = Shat.start,
                                                                              usePrior = FALSE,
                                                                              get.CIs = p$get.CIs,
                                                                              CI.method = "wald",
                                                                              run.optimx = p$run.optimx
                                ),
                                .rep.res = rep.res )
      
      
      
      cat("\n doParallel flag: Done csm-mle-sd if applicable")
    }
    
    
    
    # ~~ CSM MCMC ------------------------------
    
    if ( "csm-mcmc" %in% all.methods ) {
      # # temp for refreshing code
      # path = "/home/groups/manishad/SAPH"
      # setwd(path)
      # source("helper_SAPH.R")
      
      # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
      #  returns 2 lines of output, one for posterior mean and one for posterior median
      # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
      rep.res = run_method_safe(method.label = c("csm-mcmc-pmean",
                                                 "csm-mcmc-pmed",
                                                 "csm-mcmc-max-lp-iterate"),
                                method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dp.csm$yi,
                                                                                   .sei = sqrt(dp.csm$vi),
                                                                                   .tcrit = dp.csm$tcrit,
                                                                                   
                                                                                   .Mu.start = Mhat.start,
                                                                                   .Tt.start = max(0.01, Shat.start),
                                                                                   .stan.adapt_delta = p$stan.adapt_delta,
                                                                                   .stan.maxtreedepth = p$stan.maxtreedepth),
                                .rep.res = rep.res )
      
      cat("\n doParallel flag: Done csm-mcmc if applicable")
    }
    
    # ~ Secondary/Sanity-Check Methods ------------------------------
    # ~~ 2PSM With CSM Dataset ------------------------------
    
    if ( "2psm-csm-dataset" %in% all.methods &
         nrow(dp.csm) > 0 ) {
      
      rep.res = run_method_safe(method.label = c("2psm-csm-dataset"),
                                method.fn = function() {
                                  mod = weightfunct( effect = dp.csm$yi,
                                                     v = dp.csm$vi,
                                                     steps = c(0.025, 1),
                                                     table = TRUE ) 
                                  
                                  H = mod[[2]]$hessian
                                  ses = sqrt( diag( solve(H) ) )
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = mod[[2]]$par[2],
                                                            MLo = mod[[2]]$par[2] - qnorm(.975) * ses[2],
                                                            MHi = mod[[2]]$par[2] + qnorm(.975) * ses[2],
                                                            
                                                            Shat = sqrt( mod[[2]]$par[1] ),
                                                            # truncate lower limit at 0
                                                            SLo = sqrt( max( 0, mod[[2]]$par[1] - qnorm(.975) * ses[1] ) ),
                                                            SHi = sqrt( mod[[2]]$par[1] + qnorm(.975) * ses[1] ) ) ) 
                                },
                                .rep.res = rep.res )
      
      cat("\n doParallel flag: Done 2psm-csm-dataset if applicable")
      
    }
    
    # ~~ Naive (Unhacked Only) ------------------------------
    
    if ( "prereg-naive" %in% all.methods &
         nrow(dp.unhacked) > 0 ) {
      rep.res = run_method_safe(method.label = c("prereg-naive"),
                                method.fn = function() {
                                  mod = rma( yi = dp.unhacked$yi,
                                             vi = dp.unhacked$vi,
                                             method = "REML",
                                             knha = TRUE )
                                  
                                  report_meta(mod, .mod.type = "rma")
                                },
                                .rep.res = rep.res )
      
      cat("\n doParallel flag: Done prereg-naive if applicable")
    }
    
    
    
    # ~~ LTMA: MLE *with only* affirms (SD param) ------------------------------
    
    # include only affirmatives
    if ( "ltn-mle-sd" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("ltn-mle-sd"),
                                method.fn = function() estimate_jeffreys_RTMA(yi = dpa$yi,
                                                                              sei = sqrt(dpa$vi),
                                                                              par2is = "Tt",
                                                                              tcrit = qnorm(0.975), 
                                                                              Mu.start = Mhat.start,
                                                                              par2.start = Shat.start,
                                                                              usePrior = FALSE,
                                                                              get.CIs = p$get.CIs,
                                                                              CI.method = "wald",
                                                                              run.optimx = p$run.optimx),
                                .rep.res = rep.res )
    }
    
    # ~~ Sanity checks: Prior and NLL Agreement Between Stan and R ------------------------------
    
    if ( FALSE ) {
      
      # ~~~ Sanity Check #1A: Agreement at Fixed Params, Nonaffirms Only -----------
      # these are calculated for a FIXED (mu=2, tau=2) to match the sanity checks embedded in 
      #  stan.model above (see "generated quantities")
      # these calls are as in nlpost_jeffreys_RTMA
      log.lkl.sanity.R = -1 * joint_nll_2( .yi = dpn$yi,
                                           .sei = sqrt(dpn$vi),
                                           .tcrit = dpn$tcrit,
                                           .Mu = 2,
                                           .Tt = 2 )
      
      
      log.prior.sanity.R = lprior( .sei = sqrt(dpn$vi),
                                   .Mu = 2,
                                   .Tt = 2,
                                   .tcrit = dpn$tcrit )
      
      #wastefully re-run MCMC in order to capture its full output (which isn't preserved
      #  when it's run inside run_method_safe)
      # with only nonaffirms:
      temp = estimate_jeffreys_mcmc_RTMA(.yi = dpn$yi,
                                         .sei = sqrt(dpn$vi),
                                         .tcrit = dpn$tcrit,
                                         .Mu.start = Mhat.start,
                                         .Tt.start = Shat.start,
                                         .stan.adapt_delta = p$stan.adapt_delta,
                                         .stan.maxtreedepth = p$stan.maxtreedepth)
      ext = rstan::extract(temp$post)
      # "unique" because these have a fixed value across all iterates
      log.lkl.sanity.stan = unique(ext$log_lik_sanity)
      log.prior.sanity.stan = unique(ext$log_prior_sanity)
      
      expect_equal( as.numeric(log.lkl.sanity.R),
                    as.numeric(log.lkl.sanity.stan),
                    tolerance = 0.001 )
      expect_equal( as.numeric(log.prior.sanity.R),
                    as.numeric(log.prior.sanity.stan),
                    tolerance = 0.001 )
      
      
      ### Sanity Check #1B: Agreement at MaxLPIterate Params ###
      # look at nlpost at the maxLPiterate vs. the MAP
      best.ind = which.max(ext$log_post)
      ( log.lkl.stan = ext$log_lik[best.ind] )
      ( log.prior.stan = ext$log_prior[best.ind] )
      
      ( log.lkl.R = -1 * joint_nll_2( .yi = dpn$yi,
                                      .sei = sqrt(dpn$vi),
                                      .tcrit = dpn$tcrit,
                                      .Mu = ext$mu[best.ind],
                                      .Tt = ext$tau[best.ind] ) )
      
      ( log.prior.R = lprior( .sei = sqrt(dpn$vi),
                              .Mu = ext$mu[best.ind],
                              .Tt = ext$tau[best.ind],
                              .tcrit = dpn$tcrit ) )
      
      # THESE AGREE!
      expect_equal( as.numeric(log.lkl.R),
                    as.numeric(log.lkl.stan),
                    tolerance = 0.001 )
      expect_equal( as.numeric(log.prior.R),
                    as.numeric(log.prior.stan),
                    tolerance = 0.001 )
      
      
      ( log.post.R = -1 * nlpost_jeffreys_RTMA(.pars = c(ext$mu[best.ind],
                                                         ext$tau[best.ind]),
                                               .par2is = "Tt",
                                               .yi = dpn$yi,
                                               .sei = sqrt(dpn$vi),
                                               .tcrit = dpn$tcrit,
                                               .usePrior = TRUE) )
      
      
      ( log.post.stan = ext$log_post[best.ind] )
      
      # MATCHES :)
      expect_equal( as.numeric(log.post.R),
                    as.numeric(log.post.stan),
                    tolerance = 0.001 )
      
      
      
      
      
      # ~~~ Sanity Check #2A: Agreement at Fixed Params, CSM Dataset (Includes Affirms) -------
      # these are calculated for a FIXED (mu=2, tau=2) to match the sanity checks embedded in 
      #  stan.model above (see "generated quantities")
      # these calls are as in nlpost_jeffreys_RTMA
      
      log.lkl.sanity.R = -1 * joint_nll_2( .yi = dp.csm$yi,
                                           .sei = sqrt(dp.csm$vi),
                                           .tcrit = dp.csm$tcrit,
                                           .Mu = 2,
                                           .Tt = 2 )
      
      #wastefully re-run MCMC in order to capture its full output (which isn't preserved
      #  when it's run inside run_method_safe)
      # with only nonaffirms:
      temp = estimate_jeffreys_mcmc_RTMA(.yi = dp.csm$yi,
                                         .sei = sqrt(dp.csm$vi),
                                         .tcrit = dp.csm$tcrit,
                                         .Mu.start = Mhat.start,
                                         .Tt.start = Shat.start,
                                         .stan.adapt_delta = p$stan.adapt_delta,
                                         .stan.maxtreedepth = p$stan.maxtreedepth)
      ext = rstan::extract(temp$post)
      # "unique" because these have a fixed value across all iterates
      log.lkl.sanity.stan = unique(ext$log_lik_sanity)
      
      # agrees :)
      expect_equal( as.numeric(log.lkl.sanity.R),
                    as.numeric(log.lkl.sanity.stan),
                    tolerance = 0.001 )
      
      # not checking prior because helper_SAPH::lprior doesn't handle affirms
      
      
      
      # ~~~ Sanity Check #3A: Agreement at Fixed Params, Affirms Only --------
      
      dpa = dp %>% filter(affirm == TRUE)
      
      
      ( log.lkl.sanity.R = -1 * joint_nll_2( .yi = dpa$yi,
                                             .sei = sqrt(dpa$vi),
                                             .tcrit = dpa$tcrit,
                                             .Mu = 2,
                                             .Tt = 2 ) )
      
      
      # write my own simpler lkl
      library(truncnorm)
      log.lkl.sanity.R.2 = sum( log( dtruncnorm( x = dpa$yi,
                                                 mean = 2,
                                                 sd = sqrt(2^2 + dpa$vi),
                                                 a = dpa$tcrit * sqrt(dpa$vi),
                                                 b = 99 ) ) )
      
      expect_equal(log.lkl.sanity.R, log.lkl.sanity.R.2)
      
      # another one, closer to init_stan_model's generated_quantities block
      term1 = dnorm(x = dpa$yi,
                    mean = 2,
                    sd = sqrt(2^2 + dpa$vi),
                    log = TRUE )
      critScaled = dpa$tcrit * sqrt(dpa$vi)
      term2 = log( 1 - pnorm(q = critScaled,
                             mean = 2,
                             sd = sqrt(2^2 + dpa$vi) ) )
      # agrees!!!!!
      expect_equal( sum( term1 - term2 ),
                    log.lkl.sanity.R )
      
      
      #wastefully re-run MCMC in order to capture its full output (which isn't preserved
      #  when it's run inside run_method_safe)
      # with only nonaffirms:
      temp = estimate_jeffreys_mcmc_RTMA(.yi = dpa$yi,
                                         .sei = sqrt(dpa$vi),
                                         .tcrit = dpa$tcrit,
                                         .Mu.start = Mhat.start,
                                         .Tt.start = Shat.start,
                                         .stan.adapt_delta = p$stan.adapt_delta,
                                         .stan.maxtreedepth = p$stan.maxtreedepth)
      ext = rstan::extract(temp$post)
      # "unique" because these have a fixed value across all iterates
      ( log.lkl.sanity.stan = unique(ext$log_lik_sanity) )
      
      # agrees
      expect_equal( as.numeric(log.lkl.sanity.R),
                    as.numeric(log.lkl.sanity.stan),
                    tolerance = 0.001 )
      
      # not checking prior because helper_SAPH::lprior doesn't handle affirms
      
      # # # if needed to refresh code
      # path = "/home/groups/manishad/SAPH"
      # setwd(path)
      # source("helper_SAPH.R")
    }
    
    # # SAVE 
    # # methods from earlier simulations where I was bias-correcting the affirmatives
    # # # ~~ Bias-Corrected Estimator #1
    # # omniscient version: we know which studies are unhacked
    # #  (also includes some affirmatives)
    # 
    # corr1 = correct_dataset_phack(.dp = dp,
    #                               .p = p,
    #                               hackAssumption = "omniscient")
    # 
    # expect_equal( corr1$sanityChecks$kAssumedHacked,
    #               sum(dp$hack == p$hack) )
    # 
    # expect_equal( corr1$modUH$k,
    #               sum(dp$hack == "no") )
    # 
    # # add to results
    # repRes = add_method_result_row(repRes = NA,
    #                                corrObject = corr1,
    #                                methName = "omniscient")
    # 
    # # ~~ Bias-Corrected Estimator #2
    # # omniscient version: we know which studies are unhacked
    # #  (also includes some affirmatives)
    # 
    # corr2 = correct_dataset_phack(.dp = dp,
    #                               .p = p,
    #                               hackAssumption = "allAffirms")
    # 
    # # this check assumes p$hack = "affirm" instead of "signif"
    # expect_equal( corr2$sanityChecks$kAssumedHacked,
    #               sum(dp$affirm) )
    # 
    # expect_equal( corr2$modUH$k,
    #               sum(dp$affirm == FALSE) )
    # 
    # 
    # # add to results
    # repRes = add_method_result_row(repRes = repRes,
    #                                corrObject = corr2,
    #                                methName = "allAffirms")
    
    
    # # SAVE: OLD RETRUN STRUCTURE
    # # ~ Write Results
    # 
    # # add stats for the simulated dataset that don't change based on correction method
    # # these are prefaced by dataset name for clarity
    # # **note: all columns before methName don't depend on the correction method used
    # repRes = repRes %>% add_column(repName = i,
    #                                
    #                                # dp.k = nrow(dp),
    #                                # dp.kAffirm = sum(dp$affirm == TRUE),
    #                                # dp.kNonaffirm = sum(dp$affirm == FALSE),
    #                                # dp.Nrealized = mean(dp$N),
    #                                
    #                                report_rma(modAll,
    #                                           .Mu = p$Mu,
    #                                           .suffix = "All"),
    #                                report_rma(modPub,
    #                                           .Mu = p$Mu,
    #                                           .suffix = "Naive"),
    #                                
    #                                .before = 1 )
    # 
    # # add in scenario parameters
    # repRes = repRes %>% add_column( scenName = scen,
    #                                 p,
    #                                 .after = 1 )
    # 
    # 
    # repRes
    # 
    # cat("\nSURVIVED MAKING REPRES:")
    # print(repRes)
    
    # ~ Add Scen Params and Sanity Checks
    
    # add in scenario parameters
    # do NOT use rbind here; bind_cols accommodates possibility that some methods' rep.res
    #  have more columns than others
    rep.res = p %>% bind_cols( rep.res )
    
    # add more info
    rep.res = rep.res %>% add_column( rep.name = i, .before = 1 )
    rep.res = rep.res %>% add_column( scen.name = scen, .before = 1 )
    rep.res = rep.res %>% add_column( job.name = jobname, .before = 1 )
    
    
    cat("\ndoParallel flag: Before adding sanity checks to rep.res")
    
    
    # some san.checks will fail for sim.env = stefan b/c e.g., we don't have within-study sample sizes
    if ( p$sim.env == "mathur" ) {
      # add info about simulated datasets
      # "ustudies"/"udraws" refers to underlying studies/draws prior to hacking or publication bias
      ( sancheck.prob.ustudies.published =  mean( d.first$study %in% unique(dp$study) ) )
      expect_equal( sancheck.prob.ustudies.published, nrow(dp)/nrow(d.first) )
      # this one should always be 100% unless there's also publication bias:
      ( sancheck.prob.unhacked.ustudies.published =  mean( d.first$study[ d.first$hack == "no" ] %in% unique( dp$study[ dp$hack == "no" ] ) ) )
      # under affim hacking, will be <100%:
      ( sancheck.prob.hacked.ustudies.published =  mean( d.first$study[ d.first$hack != "no" ] %in% unique( dp$study[ dp$hack != "no" ] ) ) )
      
      # might NOT be 100% if you're generating multiple draws per unhacked studies but favoring, e.g., a random one:
      ( sancheck.prob.unhacked.udraws.published =  mean( d$study.draw[ d$hack == "no" ] %in% unique( dp$study.draw[ dp$hack == "no" ] ) ) )
      ( sancheck.prob.hacked.udraws.published =  mean( d$study.draw[ d$hack != "no" ] %in% unique( dp$study.draw[ dp$hack != "no" ] ) ) )
      
      
      
      #*this one is especially important: under worst-case hacking, it's analogous to prop.retained  in
      #  TNE since it's the proportion of the underlying distribution that's nonaffirmative
      ( sancheck.prob.unhacked.udraws.nonaffirm =  mean( d$affirm[ d$hack == "no" ] == FALSE ) )
      # a benchmark for average power:
      ( sancheck.prob.unhacked.udraws.affirm =  mean( d$affirm[ d$hack == "no" ] ) )
      ( sancheck.prob.hacked.udraws.nonaffirm =  mean( d$affirm[ d$hack != "no" ] == FALSE ) )
      ( sancheck.prob.hacked.udraws.affirm =  mean( d$affirm[ d$hack != "no" ] ) )
      
      # probability that a published, nonaffirmative draw is from a hacked study
      # under worst-case hacking, should be 0
      ( sancheck.prob.published.nonaffirm.is.hacked = mean( dp$hack[ dp$affirm == 0 ] != "no" ) )
      # this will be >0
      ( sancheck.prob.published.affirm.is.hacked = mean( dp$hack[ dp$affirm == 1 ] != "no" ) )
      
      # average yi's 
      
      rep.res = rep.res %>% add_column(   sancheck.dp.k = nrow(dp),
                                          sancheck.dp.k.affirm = sum(dp$affirm == TRUE),
                                          sancheck.dp.k.nonaffirm = sum(dp$affirm == FALSE),
                                          
                                          sancheck.dp.k.affirm.unhacked = sum(dp$affirm == TRUE & dp$hack == "no"),
                                          sancheck.dp.k.affirm.hacked = sum(dp$affirm == TRUE & dp$hack != "no"),
                                          sancheck.dp.k.nonaffirm.unhacked = sum(dp$affirm == FALSE & dp$hack == "no"),
                                          sancheck.dp.k.nonaffirm.hacked = sum(dp$affirm == FALSE & dp$hack != "no"),
                                          
                                          # means draws per HACKED, published study
                                          sancheck.dp.meanN.hacked = mean( dp$N[dp$hack != "no"] ),
                                          sancheck.dp.q90N.hacked = quantile( dp$N[dp$hack != "no"], 0.90 ),
                                          
                                          # average yi's of published draws from each study type
                                          sancheck.mean.yi.unhacked.pub.study = mean( dp$yi[ dp$hack == "no"] ),
                                          sancheck.mean.yi.hacked.pub.study = mean( dp$yi[ dp$hack != "no"] ),
                                          
                                          
                                          sancheck.mean.mui.unhacked.pub.nonaffirm = mean( dp$mui[ dp$hack == "no" & dp$affirm == FALSE ] ),
                                          sancheck.mean.yi.unhacked.pub.nonaffirm = mean( dp$yi[ dp$hack == "no" & dp$affirm == FALSE ] ),
                                          sancheck.mean.yi.unhacked.pub.affirm = mean( dp$yi[ dp$hack == "no" & dp$affirm == TRUE ] ),
                                          
                                          sancheck.mean.yi.hacked.pub.nonaffirm = mean( dp$yi[ dp$hack != "no" & dp$affirm == FALSE ] ),
                                          sancheck.mean.yi.hacked.pub.affirm = mean( dp$yi[ dp$hack != "no" & dp$affirm == TRUE ] ),
                                          
                                          # average Zi's
                                          sancheck.mean.Zi.unhacked.pub.study = mean( dp$Zi[ dp$hack == "no"] ),
                                          sancheck.mean.Zi.hacked.pub.study = mean( dp$Zi[ dp$hack != "no"] ),
                                          
                                          sancheck.mean.Zi.unhacked.pub.nonaffirm = mean( dp$Zi[ dp$hack == "no" & dp$affirm == FALSE ] ),
                                          sancheck.mean.Zi.unhacked.pub.affirm = mean( dp$Zi[ dp$hack == "no" & dp$affirm == TRUE ] ),
                                          
                                          sancheck.mean.Zi.hacked.pub.nonaffirm = mean( dp$Zi[ dp$hack != "no" & dp$affirm == FALSE ] ),
                                          sancheck.mean.Zi.hacked.pub.affirm = mean( dp$Zi[ dp$hack != "no" & dp$affirm == TRUE ] ),
                                          
                                          
                                          sancheck.prob.ustudies.published = sancheck.prob.ustudies.published,
                                          sancheck.prob.unhacked.ustudies.published = sancheck.prob.unhacked.ustudies.published,
                                          sancheck.prob.hacked.ustudies.published = sancheck.prob.hacked.ustudies.published,
                                          
                                          sancheck.prob.unhacked.udraws.published = sancheck.prob.unhacked.udraws.published,
                                          sancheck.prob.hacked.udraws.published = sancheck.prob.hacked.udraws.published,
                                          
                                          sancheck.prob.unhacked.udraws.nonaffirm = sancheck.prob.unhacked.udraws.nonaffirm,
                                          sancheck.prob.unhacked.udraws.affirm = sancheck.prob.unhacked.udraws.affirm,
                                          sancheck.prob.hacked.udraws.nonaffirm = sancheck.prob.hacked.udraws.nonaffirm,
                                          sancheck.prob.hacked.udraws.affirm = sancheck.prob.hacked.udraws.affirm,
                                          
                                          sancheck.prob.published.nonaffirm.is.hacked = sancheck.prob.published.nonaffirm.is.hacked
      )
    }
    
    rep.res
    
  }  ### end foreach loop
  
} )[3]  # end system.time


# quick look
rs %>% dplyr::select(method, Mhat, MLo, MHi)




table(rs$method)


# # LOCAL
# # how badly biased are the hacked studies?
# temp = rs %>% filter(method == "naive")
# mean(temp$sancheck.mean.yi.hacked.pub.study)
# mean(temp$sancheck.mean.yi.unhacked.pub.study)
# 
# mean(temp$sancheck.mean.yi.unhacked.pub.affirm)
# mean(temp$sancheck.mean.yi.hacked.pub.affirm)
# 
# mean(temp$sancheck.mean.yi.unhacked.pub.nonaffirm)
# mean(temp$sancheck.mean.yi.hacked.pub.nonaffirm)



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

#@debugging
#rs$overall.error[rs$method == "rtma-pkg"]


# ~ QUICK RESULTS SUMMARY ---------------------------------------------------

if ( run.local == TRUE ) {
  
  # quick look locally
  # rs %>% mutate_if(is.numeric, function(x) round(x,2) )
  
  agg = rs %>% group_by(method) %>%
    summarise( PropMhatNA = mean(is.na(Mhat)),
               PropCI.NA = mean(is.na(MLo)),
               
               Mhat = meanNA(Mhat),
               MhatMSE = meanNA( (Mhat - Mu)^2 ),
               MhatBias = meanNA(Mhat - Mu),
               MhatEmpSE = sd( Mhat, na.rm = TRUE ),
               #ShatMn = meanNA(Shat),
               
               MhatCover = meanNA(MLo < Mu & MHi > Mu),
               MhatWidth = meanNA( MHi - MLo ),
               MLo = meanNA(MLo),
               MHi = meanNA(MHi) )
  
  # round
  agg = as.data.frame( agg %>% mutate_if( is.numeric,
                                          function(x) round(x,2) ) )
  
  agg
  
  
  # scenario diagnostics for scenario
  keepers = namesWith("sancheck.", rs)
  agg.checks = rs %>% summarise_at( keepers,
                                    function(x) round( mean(x), 2) )
  
  t(agg.checks)
  
}




# ~ WRITE LONG RESULTS ------------------------------
if ( run.local == FALSE ) {
  setwd("/home/groups/manishad/SAPH/long_results")
  fwrite( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
}