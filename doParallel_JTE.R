
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
# path = "/home/groups/manishad/JTE"
# setwd(path)
# source("doParallel_JTE.R")


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
  source("stefan_phackR_fns.R")
  
  
  # ~~ ****** Set Local Sim Params -----------------------------
  
  scen.params = tidyr::expand_grid(
    # full list (save):
    rep.methods = "REML ; ML ; DL ; PMM ; EB ; robu ; jeffreys",
    
    # *If you reorder the args, need to adjust wrangle_agg_local
    ### args shared between sim environments
    k.pub = c(50),  # intentionally out of order so that jobs with boundary choices with complete first
    hack = c("affirm"),
    prob.hacked = c(0),
    # important: if sim.env = stefan, these t2 args are ONLY used for setting start values
    #   and for checking bias of Shat, so set them to have the correct t2a
    #   not clear what t2w should be given the way stefan implements hacking
    t2a = c(0.2^2),
    t2w = c(0),
    # same with Mu
    Mu = c(0.5),
    
    Nmax = 1,
    m = 50,
    true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
    rho = c(0),
    
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
  
  setwd("/home/groups/manishad/JTE/applied_examples/data")
  d.lodder = fread("lodder_prepped.csv")
  lodder.ses = sqrt(d.lodder$vi)
  
  cat("\n\ndoParallel: just read in Lodder SEs:")
  summary(lodder.ses)
  
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
    
    
    d = sim_meta_2( Nmax = p$Nmax,
                    Mu = p$Mu,
                    t2a = p$t2a,
                    m = p$m,
                    t2w = p$t2w,
                    true.sei.expr = p$true.sei.expr,
                    hack = p$hack,
                    rho = p$rho,
                    
                    k.pub = p$k.pub,
                    prob.hacked = p$prob.hacked,
                    return.only.published = FALSE)
    
    
    
    d$Zi = d$yi / sqrt(d$vi)
    d$sei = sqrt(d$vi)
    
    if ( i == 1 ) cat("\n\nHEAD OF D:\n")
    if ( i == 1 ) print(head(d))
    
    
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
    
    # ~~ Metafor heterogeneity estimators ------------------------------
  
    # pg 282:
    # https://cran.r-project.org/web/packages/metafor/metafor.pdf
    metafor.methods = all.methods[ all.methods %in% c("REML", "ML", "DL", "EB", "PMM", "HS", "SJ") ]
    

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
    
    
    # ~~ Robust Variance Estimation ------------------------------
    
    
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
    # ~~ ***** MCMC ------------------------------
    
    if ( "jeffreys" %in% all.methods ) {
      # # temp for refreshing code
      # path = "/home/groups/manishad/JTE"
      # setwd(path)
      # source("helper_JTE.R")
      # source("init_stan_model_JTE.R")
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
      
      
      Mhat.MaxLP = rep.res$Mhat[ rep.res$method == "jeffreys-max-lp-iterate" ]
      Shat.MaxLP = rep.res$Shat[ rep.res$method == "jeffreys-max-lp-iterate" ]
      
      cat("\n doParallel flag: Done jeffreys-mcmc if applicable")
    }
    
    # NOTE: if doing run.local, this will break if you didn't run naive
    if (run.local == TRUE) srr(rep.res)
    
    
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
    if (FALSE) {
      # add info about simulated datasets
      # "ustudies"/"udraws" refers to underlying studies/draws prior to hacking or publication bias
      ( sancheck.prob.ustudies.published =  mean( d.first$study %in% unique(d$study) ) )
      expect_equal( sancheck.prob.ustudies.published, nrow(dp)/nrow(d.first) )
      # this one should always be 100% unless there's also publication bias:
      ( sancheck.prob.unhacked.ustudies.published =  mean( d.first$study[ d.first$hack == "no" ] %in% unique( d$study[ d$hack == "no" ] ) ) )
      # under affim hacking, will be <100%:
      ( sancheck.prob.hacked.ustudies.published =  mean( d.first$study[ d.first$hack != "no" ] %in% unique( d$study[ d$hack != "no" ] ) ) )
      
      # might NOT be 100% if you're generating multiple draws per unhacked studies but favoring, e.g., a random one:
      ( sancheck.prob.unhacked.udraws.published =  mean( d$study.draw[ d$hack == "no" ] %in% unique( d$study.draw[ d$hack == "no" ] ) ) )
      ( sancheck.prob.hacked.udraws.published =  mean( d$study.draw[ d$hack != "no" ] %in% unique( d$study.draw[ d$hack != "no" ] ) ) )
      
      
      
      #*this one is especially important: under worst-case hacking, it's analogous to prop.retained  in
      #  TNE since it's the proportion of the underlying distribution that's nonaffirmative
      ( sancheck.prob.unhacked.udraws.nonaffirm =  mean( d$affirm[ d$hack == "no" ] == FALSE ) )
      # a benchmark for average power:
      ( sancheck.prob.unhacked.udraws.affirm =  mean( d$affirm[ d$hack == "no" ] ) )
      ( sancheck.prob.hacked.udraws.nonaffirm =  mean( d$affirm[ d$hack != "no" ] == FALSE ) )
      ( sancheck.prob.hacked.udraws.affirm =  mean( d$affirm[ d$hack != "no" ] ) )
      
      # probability that a published, nonaffirmative draw is from a hacked study
      # under worst-case hacking, should be 0
      ( sancheck.prob.published.nonaffirm.is.hacked = mean( d$hack[ d$affirm == 0 ] != "no" ) )
      # this will be >0
      ( sancheck.prob.published.affirm.is.hacked = mean( d$hack[ d$affirm == 1 ] != "no" ) )
      
      # average yi's 
      
      rep.res = rep.res %>% add_column(   sancheck.dp.k = nrow(dp),
                                          sancheck.dp.k.affirm = sum(d$affirm == TRUE),
                                          sancheck.dp.k.nonaffirm = sum(d$affirm == FALSE),
                                          
                                          sancheck.dp.k.affirm.unhacked = sum(d$affirm == TRUE & d$hack == "no"),
                                          sancheck.dp.k.affirm.hacked = sum(d$affirm == TRUE & d$hack != "no"),
                                          sancheck.dp.k.nonaffirm.unhacked = sum(d$affirm == FALSE & d$hack == "no"),
                                          sancheck.dp.k.nonaffirm.hacked = sum(d$affirm == FALSE & d$hack != "no"),
                                          
                                          # means draws per HACKED, published study
                                          sancheck.dp.meanN.hacked = mean( d$N[d$hack != "no"] ),
                                          sancheck.dp.q90N.hacked = quantile( d$N[d$hack != "no"], 0.90 ),
                                          
                                          # average yi's of published draws from each study type
                                          sancheck.mean.yi.unhacked.pub.study = mean( d$yi[ d$hack == "no"] ),
                                          sancheck.mean.yi.hacked.pub.study = mean( d$yi[ d$hack != "no"] ),
                                          
                                          
                                          sancheck.mean.mui.unhacked.pub.nonaffirm = mean( d$mui[ d$hack == "no" & d$affirm == FALSE ] ),
                                          sancheck.mean.yi.unhacked.pub.nonaffirm = mean( d$yi[ d$hack == "no" & d$affirm == FALSE ] ),
                                          sancheck.mean.yi.unhacked.pub.affirm = mean( d$yi[ d$hack == "no" & d$affirm == TRUE ] ),
                                          
                                          sancheck.mean.yi.hacked.pub.nonaffirm = mean( d$yi[ d$hack != "no" & d$affirm == FALSE ] ),
                                          sancheck.mean.yi.hacked.pub.affirm = mean( d$yi[ d$hack != "no" & d$affirm == TRUE ] ),
                                          
                                          # average Zi's
                                          sancheck.mean.Zi.unhacked.pub.study = mean( d$Zi[ d$hack == "no"] ),
                                          sancheck.mean.Zi.hacked.pub.study = mean( d$Zi[ d$hack != "no"] ),
                                          
                                          sancheck.mean.Zi.unhacked.pub.nonaffirm = mean( d$Zi[ d$hack == "no" & d$affirm == FALSE ] ),
                                          sancheck.mean.Zi.unhacked.pub.affirm = mean( d$Zi[ d$hack == "no" & d$affirm == TRUE ] ),
                                          
                                          sancheck.mean.Zi.hacked.pub.nonaffirm = mean( d$Zi[ d$hack != "no" & d$affirm == FALSE ] ),
                                          sancheck.mean.Zi.hacked.pub.affirm = mean( d$Zi[ d$hack != "no" & d$affirm == TRUE ] ),
                                          
                                          
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
#rs %>% dplyr::select(method, Shat, SLo, SHi, Mhat, MLo, MHi)




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
               ShatEmpSE = sd( Shat, na.rm = TRUE ),
               #ShatMn = meanNA(Shat),
               
               ShatCover = meanNA(SLo < tau & SHi > tau),
               ShatWidth = meanNA( SHi - SLo ),
               SLo = meanNA(SLo),
               SHi = meanNA(SHi),
               
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
  setwd("/home/groups/manishad/JTE/long_results")
  fwrite( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
}