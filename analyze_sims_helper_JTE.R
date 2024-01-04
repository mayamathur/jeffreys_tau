
# ~ FNS FOR SUMMARIZING SIMULATION DATA ----------------------------------------------


# taken from TNE 2022-2-27
# fn for aggregating so we can look at different
#  iterate-level filtering rules
# .s: the iterate-level stitched data (not yet aggregated in any way)
# averagefn: fn to use when aggregating results across scenarios
# expected.sim.reps: only used for sanity checks
make_agg_data = function( .s,
                          .averagefn = "median",
                          badCoverageCutoff = 0.85,
                          expected.sim.reps = NA ){
  
  # # TEST ONLY
  # .s = s
  # .averagefn = "median",
  # badCoverageCutoff = 0.85,
  # expected.sim.reps = NA
  
  
  # make unique scenario variable, defined as scen.name AND method
  if ( !"unique.scen" %in% names(.s) ) .s$unique.scen = paste(.s$sim.env, .s$scen.name, .s$method)
  
  ##### Outcome and Parameter Variables #####
  # "outcome" variables used in analysis
  analysis.vars = c( 
    "Mhat",
    "Shat",
    
    "MLo",
    "MHi",
    
    "SLo",
    "SHi",
    
    
    ##### variables to be created in mutate below:
    
    "MhatBias",
    "ShatBias",
    
    # "MhatRelBias",
    # "VhatRelBias",
    
    "MhatCover",
    "ShatCover",
    
    "MhatWidth",
    "ShatWidth",
    
    "MhatRMSE",
    "ShatRMSE",
    
    "MhatEstSE",
    "ShatEstSE",
    
    "MhatEmpSE",
    "ShatEmpSE",
    
    # diagnostics regarding point estimation and CIs
    "MhatEstFail",
    "MhatCIFail",
    "ShatEstFail",
    "ShatCIFail",
    
    # #@2022-3-10 TEMP: COMMENTED OUT BECAUSE I DIDN'T RUN OPTIMX METHODS, SO THIS BREAKS
    # # diagnostics for main Mhat optimizer
    # "OptimConverged",
    # 
    # #@2022-3-10 TEMP: COMMENTED OUT BECAUSE I DIDN'T RUN OPTIMX METHODS, SO THIS BREAKS
    # # diagnostics for other optimizers
    # "OptimxMhatWinner",
    # "OptimxPropAgreeMhatWinner",
    # "OptimxPropAgreeConvergersMhatWinner",
    # 
    # "OptimxShatWinner",
    # "OptimxPropAgreeShatWinner",
    # "OptimxPropAgreeConvergersShatWinner",
    
    # Stan diagnostics, part 2
    "StanWarned",
    "MhatRhatGt1.01",
    "MhatRhatGt1.05",
    "MhatRhatGt1.10",
    "MhatRhatMax",
    
    "ShatRhatGt1.01",
    "ShatRhatGt1.05",
    "ShatRhatGt1.10",
    "ShatRhatMax"
  )
  
  
  
  
  # variables that define the scenarios
  
  
  if ( all(.s$sim.env == "mathur") ) {
    param.vars = c("unique.scen",  
                   "method",
                   "sim.env",
                   
                   "Nmax",
                   "Mu",
                   "t2a",
                   "t2w",
                   "m",
                   "true.sei.expr",
                   "hack",
                   "rho",
                   "k.pub.nonaffirm",
                   "prob.hacked",
                   "stan.adapt_delta",
                   "stan.maxtreedepth")
    
  } else if ( all(.s$sim.env == "stefan") ) {
    param.vars = c("unique.scen",  
                   "method",
                   "sim.env",
                   
                   "Mu",
                   "t2a",
                   "t2w",
                   
                   "hack",
                   "strategy.stefan",
                   "alternative.stefan",
                   "stringent.hack",
                   
                   "k.pub.nonaffirm",
                   "prob.hacked",
                   "stan.adapt_delta",
                   "stan.maxtreedepth")
  }
  
 # TEMP aggregation issue
  param.vars = c("unique.scen",  
                 "method",
                 "sim.env",
                 
                 "Nmax",
                 "Mu",
                 "t2a",
                 "t2w",
                 "m",
                 "true.sei.expr",
                 "hack",
                 "rho",
                 "k.pub.nonaffirm",
                 "prob.hacked",
                 "stan.adapt_delta",
                 "stan.maxtreedepth",
                 
                 "strategy.stefan",
                 "alternative.stefan",
                 "stringent.hack")
  
  
  # sanity check to make sure we've listed all param vars
  t = .s %>% group_by_at(param.vars) %>% summarise(n())
  if ( !is.na(expected.sim.reps) ) {
    if ( max(t$`n()`) > expected.sim.reps ) stop("param.vars in make_agg_data might be missing something because grouping that way indicated some scenarios had more than expected.sim.reps")
  }
  
  
  ##### Overwrite Analysis Variables As Their Within-Scenario Means #####
  # organize variables into 3 mutually exclusive sets: 
  # - param.vars: parameter variables for grouping
  # - toDrop: variables to drop completely
  # - firstOnly: variables that are static within a scenario, for which we
  #   should just take the first one (but not param variables)
  # - takeMean: variables for which we should take the mean within scenarios
  
  #names(.s)[ !names(.s) %in% param.vars ]  # look at names of vars that need categorizing
  
  .s$V = .s$t2a + .s$t2w
  .s$S = sqrt(.s$t2a + .s$t2w)
  

  toDrop = c("rep.methods",
             "get.CIs",
             "error",
             "rep.name",
             #"doParallel.seconds",
             "optim.converged",
             "stan.warned",
             "stan.warning", #@ADDED FOR STEFAN
             "job.name",
             "overall.error",  #@ADDED FOR STEFAN
             names_with(.dat = .s, .pattern = "optimx") )
  
  firstOnly = c("scen.name",
                "unique.scen",
                "V",  # calculated from scen params
                "S")
  
  ##### Add New Variables Calculated at Scenario Level #####
  
  # if you have 10K iterates, script breaks from here forward if running locally
  # "vector memory limits"
  s2 = .s %>%
    rename(
      # static within scenario
      # just renaming for clarity
      MhatEstSE = MhatSE,
      ShatEstSE = ShatSE ) %>%
    
    # take just first entry of non-parameter variables that are static within scenarios
    group_by_at(param.vars) %>%
    mutate_at( firstOnly, 
               function(x) x[1] ) %>%
    
    # make certain ad hoc variables that don't conform to below rules
    # this step creates variables that are repeated for every rep within 
    #  a scenario, which is intentional
    
    # make variables that are calculated within scenarios
    # some of the vars are defined at the iterate level (i.e., they still vary within scen), 
    #  while others are calculated at the scen level (i.e., they are static within scen)
    # after this step, we take the means within scens of all these vars, which is immaterial
    #   for the ones that are already static within scenario
    group_by_at(param.vars) %>%
    
    mutate( sim.reps.actual = n(),
            
            # varies within scenario
            MhatBias = Mhat - Mu,
            MhatAbsBias = abs(Mhat - Mu),
            ShatBias = Shat - S,
            
            # varies within scenario
            MhatCover = covers(truth = Mu, lo = MLo, hi = MHi),
            ShatCover = covers(truth = S, lo = SLo, hi = SHi),
            
            # varies within scenario
            MhatWidth = MHi - MLo,
            ShatWidth = SHi - SLo,
            
            # varies within scenario
            MhatTestReject = MLo > 0,
            
            # static within scenario
            MhatRMSE = sqrt( meanNA( (Mhat - Mu)^2 ) ),
            ShatRMSE = sqrt( meanNA( ( Shat - S )^2 ) ),
            
            # static within scenario
            MhatEstFail = mean(is.na(Mhat)),
            MhatCIFail = mean(is.na(MLo)),
            ShatEstFail = mean(is.na(Shat)),
            ShatCIFail = mean(is.na(SLo)),
            
            # static within scenario
            MhatEmpSE = sd(Mhat, na.rm = TRUE),
            ShatEmpSE = sd(Shat, na.rm = TRUE),
            
            
            # varies within scenario
            # how much smaller is estimated SE compared to empirical one?
            MhatSEBias = MhatEstSE - MhatEmpSE,
            ShatSEBias = ShatEstSE - ShatEmpSE,
            
            # varies within scenario
            MhatSERelBias = (MhatEstSE - MhatEmpSE) / MhatEmpSE, 
            ShatSERelBias = (ShatEstSE - ShatEmpSE) / ShatEmpSE,
            
            # static within scenario
            # # ALSO COMMENTED OUT PART OF ANALYSIS.VARS ABOVE
            # OptimConverged = meanNA(optim.converged),
            # OptimxNConvergers = meanNA(optimx.Nconvergers),
            # OptimxNAgreeOfConvergersMhatWinner = meanNA(optimx.Nagree.of.convergers.Mhat.winner),
            # 
            # OptimxMhatWinner = meanNA(optimx.Mhat.winner),
            # OptimxPropAgreeMhatWinner = meanNA(optimx.Pagree.Mhat.winner),
            # OptimxPropAgreeConvergersMhatWinner = meanNA(optimx.Pagree.of.convergers.Mhat.winner),
            # 
            # OptimxShatWinner = meanNA(optimx.Shat.winner),
            # OptimxPropAgreeShatWinner = meanNA(optimx.Pagree.Shat.winner),
            # OptimxPropAgreeConvergersShatWinner = meanNA(optimx.Pagree.of.convergers.Shat.winner),
            
            # static within scenario
            StanWarned = meanNA(stan.warned),
            MhatRhatGt1.01 = meanNA(MhatRhat > 1.01),
            MhatRhatGt1.05 = meanNA(MhatRhat > 1.05),
            MhatRhatGt1.10 = meanNA(MhatRhat > 1.10),
            MhatRhatMax = max(MhatRhat),
            
            ShatRhatGt1.01 = meanNA(ShatRhat > 1.01),
            ShatRhatGt1.05 = meanNA(ShatRhat > 1.05),
            ShatRhatGt1.10 = meanNA(ShatRhat > 1.10),
            ShatRhatMax = max(ShatRhat),
            
            # SLURM timing stats
            doParallelSeconds = meanNA(doParallel.seconds),
            # minor note: even within scens, doParallel.seconds is repeated
            #  for every sim rep within the scen
            doParallelSecondsQ95 = quantile(doParallel.seconds,
                                            0.95, na.rm = TRUE),
    ) 
  
  
  # now look for which variables should have their means taken
  # this step must happen here, after we've started making s2, 
  #  so that the takeMean vars are actually in s2
  ( takeMean = names(s2)[ !names(s2) %in% c(param.vars, toDrop, firstOnly) ] )
  # sanity check: have all variables been sorted into these categories?
  expect_equal( TRUE,
                all( names(s2) %in% c(param.vars, toDrop, firstOnly, takeMean) ) )
  
  
  ##### Aggregate to Scenario Level #####
  
  # calculate scenario-level averages, but keep dataset at the rep level
  #  for now to facilitate sanity checks
  # IMPORTANT: this uses meanNA regardless of the passed .avgfun 
  #  because right now we are only aggregating WITHIN scens
  #  so we should never use median 
  
  # don't try to drop vars that don't exist
  toDrop = toDrop[ toDrop %in% names(s2) ]
  
  
  
  s3 = s2 %>%
    # take averages of numeric variables
    group_by_at(param.vars) %>%
    mutate_at( takeMean,
               function(x) meanNA(x) ) %>%
    select( -all_of(toDrop) )
  
  
  
  # sanity check: name mismatches
  trouble.vars = analysis.vars[ !analysis.vars %in% names(s2) ]
  if ( length( trouble.vars ) > 0 ) {
    stop( paste("Might have name mismatches; edit analysis.vars in make_agg_data; trouble vars are: ", trouble.vars ) )
  }
  
  # sanity check: SDs of all analysis variables should be 0 within unique scenarios
  t = data.frame( s3 %>% group_by(unique.scen) %>%
                    summarise_at( analysis.vars, sd ) )
  
  t = t %>% select(-unique.scen)
  expect_equal( FALSE,
                any( !as.matrix( t[, 2:(ncol(t)) ] ) %in% c(0, NA, NaN) ) )
  # end sanity checks
  
  
  ##### Aggregate Data at Scenario Level #####
  # make aggregated data by keeping only first row for each 
  #  combination of scenario name and calib.method
  agg = s3[ !duplicated(s3$unique.scen), ]
  
  ##### Create Variables That Are Defined At Scenario Rather Than Iterate Level #####
  agg = agg %>% mutate( BadMhatCover = MhatCover < badCoverageCutoff,
                        BadShatCover = ShatCover < badCoverageCutoff )
  
  return(agg %>% ungroup() )
}


# This is meant to be called after make_agg_data
# Can be run locally even when agg is huge
# This fn is separate from make_agg_data because it needs more frequent modification
wrangle_agg_local = function(agg) {
  ##### Make New Variables At Scenario Level ##### 
  
  # label methods more intelligently for use in plots
  agg$method.pretty = NA
  agg$method.pretty[ agg$method == c("naive") ] = "Uncorrected"
  agg$method.pretty[ agg$method == c("gold-std") ] = "Gold standard"
  agg$method.pretty[ agg$method == c("maon") ] = "MAN"
  agg$method.pretty[ agg$method == c("2psm") ] = "SM"
  agg$method.pretty[ agg$method == c("2psm-csm-dataset") ] = "SMKH" # "known hacking"
  agg$method.pretty[ agg$method == c("prereg-naive") ] = "Unhacked only"
  agg$method.pretty[ agg$method == c("pet-peese") ] = "PET-PEESE"
  agg$method.pretty[ agg$method == c("robma") ] = "RoBMA"
  #agg$method.pretty[ agg$method %in% c("jeffreys-mcmc-max-lp-iterate") ] = "RTMA"
  agg$method.pretty[ agg$method %in% c("rtma-pkg") ] = "RTMA"
  table(agg$method, agg$method.pretty)
  
  ##### Specific to Sim Env #####

  if ( "true.sei.expr" %in% names(agg) ){
    agg$true.sei.expr = as.factor(agg$true.sei.expr)
    
    agg$true.sei.expr.pretty = dplyr::recode( agg$true.sei.expr,
                                              `0.1 + rexp(n = 1, rate = 1.5)` = "sei ~ Exp(1.5)",
                                              `runif(n = 1, min = 0.1, max = 1)` = "sei ~ U(0.1, 1)",
                                              `runif(n = 1, min = 0.50, max = 0.60)` = "sei ~ U(0.5, 0.6)",
                                              `runif(n = 1, min = 0.51, max = 1.5)` = "sei ~ U(0.51, 1.5)",
                                              `runif(n = 1, min = 0.1, max = 3)` = "sei ~ U(0.1, 3)",
                                              `runif(n = 1, min = 1, max = 3)` = "sei ~ U(1, 3)",
                                              `rbeta(n = 1, 2, 5)` = "sei ~ Beta(2, 5)",
                                              `0.02 + rexp(n = 1, rate = 3)` = "sei ~ Exp(3) + 0.02",
                                              `draw_lodder_se()` = "sei from Lodder",
                                              
                                              # by default, retain original factor level
                                              .default = levels(agg$true.sei.expr) )
    print( table(agg$true.sei.expr, agg$true.sei.expr.pretty ) )
  }
  
  if ( "rho.pretty" %in% names(agg) )   agg$rho.pretty = paste("rho = ", agg$rho, sep = "")
  
  agg$MhatEstConverge = 1 - agg$MhatEstFail
  
  # indicates that RTMA is incorrectly specified
  agg$rtma.misspec = 1
  agg$rtma.misspec[ agg$hack %in% c("favor-best-affirm-wch",
                                      "affirm")] = 0
  
  return(agg)
}





# RESULTS TABLES FNS -------------------------------------------------------------

make_winner_table_col = function(.agg,
                                 yName,
                                 methods = c("naive",
                                             "gold-std",
                                             "maon",
                                             "2psm",
                                             "pcurve",
                                             "pet-peese", 
                                             "robma",
                                             "jeffreys-mcmc-pmean",
                                             "jeffreys-mcmc-pmed",
                                             "jeffreys-mcmc-max-lp-iterate", 
                                             "rtma-pkg",
                                             "prereg-naive"),
                                 summarise.fun.name = "median",
                                 digits = 2) {
  
  # # test only
  # .agg = agg
  # yName = "MhatCover"
  # methods = c("naive",
  #             "gold-std",
  #             "maon",
  #             "2psm",
  #             "pcurve",
  #             "pet-peese",
  #             "robma",
  #             "jeffreys-mcmc-pmean",
  #             "jeffreys-mcmc-pmed",
  #             "jeffreys-mcmc-max-lp-iterate",
  #             "prereg-naive")
  # summarise.fun.name = "worst10th"
  # digits = 2
  
  
  # sanity check
  if ( any( is.na( .agg$method.pretty ) ) ) {
    stop(".agg has NAs in method.pretty; will mess up group_by")
  }
  
  
  higherBetterYNames = c("MhatEstConverge")
  
  lowerBetterYNames = c("MhatAbsBias", "MhatRMSE", "MhatWidth")
  
  # Y_disp is what will be DISPLAYED in the table (e.g., retaining signs)
  .agg$Y_disp = .agg[[yName]]
  
  ##### Summarize Y_disp #####
  if ( summarise.fun.name == "mean" ) {
    .t = .agg %>% filter(method %in% methods) %>%
      group_by(method.pretty) %>%
      summarise( Y_disp = round( mean(Y_disp), digits = digits ) )
  }
  
  if ( summarise.fun.name == "median" ) {
    .t = .agg %>% filter(method %in% methods) %>%
      group_by(method.pretty) %>%
      summarise( Y_disp = round( median(Y_disp), digits = digits ) )
  }
  
  if ( summarise.fun.name == "worst10th" & yName %in% higherBetterYNames ) {
    .t = .agg %>% filter(method %in% methods) %>%
      group_by(method.pretty) %>%
      summarise( Y_disp = round( quantile(Y_disp, probs = 0.10), digits = digits ) )
  }
  
  if ( summarise.fun.name == "worst10th" & yName %in% lowerBetterYNames ) {
    .t = .agg %>% filter(method %in% methods) %>%
      group_by(method.pretty) %>%
      summarise( Y_disp = round( quantile(Y_disp, probs = 0.90), digits = digits ) )
  }
  
  # for bias, worst 10% performance is in terms of absolute value of bias
  # note this isn't abs bias
  if ( summarise.fun.name == "worst10th" & yName %in% c("MhatBias") ) {
    .t = .agg %>% filter(method %in% methods) %>%
      group_by(method.pretty) %>%
      summarise( Y_disp = round( bias_worst10th(Y_disp), digits = digits ) )
  }
  # sanity check
  #quantile(.agg$MhatBias[.agg$method.pretty == "MAN"], probs = c(0.1, .9))
  
  # MhatCover gets summarized simply as 10th percentile
  if ( summarise.fun.name == "worst10th" & yName %in% c("MhatCover") ) {
    .t = .agg %>% filter(method %in% methods) %>%
      group_by(method.pretty) %>%
      summarise( Y_disp = round( quantile(Y_disp, probs = 0.10, na.rm = TRUE), digits = digits ) )
  }
  # sanity check
  #quantile(.agg$MhatCover[.agg$method.pretty == "MAN"], probs = c(0.1))
  #quantile(.agg$MhatCover[.agg$method.pretty == "RoBMA"], probs = c(0.1))
  

  ##### Sort Best to Worst #####
  
  # Y_sort is the version that is used to compare and sort performances 
  # designed so that higher values are always better
  #  e.g., -abs(coverage-0.95)
  if ( yName %in% higherBetterYNames ) {
    .t$Y_sort = .t$Y_disp
  }
  if ( yName %in% lowerBetterYNames ) {
    .t$Y_sort = -.t$Y_disp
  }
  if ( yName %in% "MhatBias" ) {
    .t$Y_sort = -abs(.t$Y_disp)
  }
  if ( yName %in% "MhatCover" ) {
    .t$Y_sort = -abs(.t$Y_disp - 0.95)
  }
  

  # sort best to worst
  .t = .t %>% arrange( desc(Y_sort) )
  
  # remove unneeded col
  .t = .t %>% select(-Y_sort)

  
  names(.t) = c(yName, summarise.fun.name)
  .t
}

# example
# make_winner_table_col(.agg = agg,
#                       yName = "MhatCover",
#                       summarise.fun.name = "median")
# make_winner_table_col(.agg = agg,
#                       yName = "MhatCover",
#                       summarise.fun.name = "worst10th")


# for summarizing bias
bias_worst10th = function(x) {
  q10 = quantile(x, 0.10)
  q90 = quantile(x, 0.90)
  if ( abs(q10) >= abs(q90) ) return(q10)
  if ( abs(q10) < abs(q90) ) return(q90)
}
#bias_worst10th(x = rnorm(100))


# display: "xtable" or "dataframe"
make_winner_table = function( .agg,
                              .yNames = c("MhatBias",
                                          "MhatAbsBias",
                                          "MhatRMSE", 
                                          "MhatCover",
                                          "MhatWidth" ),
                                          #"MhatEstConverge"),
                              summarise.fun.name,
                              display = "dataframe"
                              #display = "xtable"
                              ){

  
  # sanity check
  if ( !( all( .yNames %in% names(.agg) ) )  ){
    not_here = paste( .yNames[ !( .yNames %in% names(.agg) ) ], collapse = ", " )
    stop( paste( "\nThe following .yNames aren't in .agg: ", not_here) )
  } 
  
  for ( .yName in .yNames ){
    newCol = make_winner_table_col(.agg = .agg,
                                   yName = .yName,
                                   summarise.fun.name = summarise.fun.name )
    
    if ( .yName == .yNames[1] ) t.all = newCol else t.all = suppressMessages( bind_cols(t.all, newCol) )
  }
  
  
  cat( paste("\n\n**** WINNER TABLE", summarise.fun.name) )
  
  cat( paste("\n\n     Number of scens:", nuni(.agg$scen.name2),
             "; proportion of all scens: ",
             round( nuni(.agg$scen.name) / nuni(agg$scen.name2), 3 ) ) )
  
  cat("\n\n")
  
  
  
  if (display == "xtable") {
    print( xtable( data.frame(t.all) ), include.rownames = FALSE )
  }
  
  if (display == "dataframe") {
    print( data.frame(t.all) )
  }
  
  #return(t.all)
  
}

# makes both winner tables (medians and worst 10th pctiles)
make_both_winner_tables = function( .agg,
                                    .yNames = c("MhatBias", "MhatAbsBias", "MhatRMSE", "MhatCover", "MhatWidth") ){
  
  make_winner_table( .agg = .agg,
                     .yNames = .yNames,
                     summarise.fun.name = "median")
  
  make_winner_table( .agg = .agg,
                     .yNames = .yNames,
                     summarise.fun.name = "mean")
  
  make_winner_table( .agg = .agg,
                     .yNames = .yNames,
                     summarise.fun.name = "worst10th")
  
}




# PLOTTING FNS -------------------------------------------------------------

# make a plot with 3 variables: x-axis, y-axis, facets, and colors
# facet vars allowed be null
quick_5var_agg_plot = function(.Xname,
                               .Yname,
                               .colorVarName,
                               .facetVar1Name = NULL,
                               .facetVar2Name = NULL,
                               
                               .dat,
                               .ggtitle = "",
                               
                               .y.breaks = NULL,
                               
                               .writePlot = FALSE,
                               .results.dir = NULL) {
  
  
  # TEST
  # agg$facetVar = paste( "rho=", agg$rho, "; ", agg$true.sei.expr.pretty, sep="")
  # table(agg$facetVar)
  # agg$rho.pretty = paste("rho = ", agg$rho, sep = "")
  # 
  # .Xname = "k.pub.nonaffirm"
  # .Yname = "MhatBias"
  # .colorVarName = "method"
  # .facetVar1Name = "rho.pretty"
  # .facetVar2Name = "true.sei.expr.pretty"
  # .dat = agg
  # .ggtitle = ""
  # .writePlot = FALSE
  # #.results.dir
  
  
  .dat$Y = .dat[[.Yname]]
  .dat$X = .dat[[.Xname]]
  .dat$colorVar = .dat[[.colorVarName]]
  # don't try to move these inside conditional statement below
  #  about facet_wrap b/c then .dat won't contain the facet vars
  .dat$facetVar1 = .dat[[.facetVar1Name]]
  .dat$facetVar2 = .dat[[.facetVar2Name]]
  
  # ~ Make base plot ----------
  p = ggplot( data = .dat,
              aes( x = X,
                   y = Y,
                   color = as.factor(colorVar) ) ) +
    
    geom_point() +
    geom_line() +
    
    # use all values of
    #scale_x_log10( breaks = unique(.dp$n) )
    # use only some values
    #scale_x_log10( breaks = c(500, 1000) ) +
    
    xlab(.Xname) +
    ylab(.Yname) +
    guides( color = guide_legend(title = .colorVarName) ) +
    ggtitle(.ggtitle) +
    theme_bw() 
  
  # ~ Add reference lines ----------
  if ( str_contains(x = .Yname, pattern = "Cover") ) {
    p = p + geom_hline( yintercept = 0.95,
                        lty = 2,
                        color = "black" ) 
    
  }
  
  if ( str_contains(x = .Yname, pattern = "Bias") ) {
    p = p + geom_hline( yintercept = 0,
                        lty = 2,
                        color = "black" ) 
    
  }
  
  # ~ Add facetting ----------
  # this block needs to be after adding geom_hlines so that the lines obey the facetting
  if ( !is.null(.facetVar1Name) & !is.null(.facetVar2Name) ) {
    p = p + facet_wrap(facetVar1 ~ facetVar2,
                       nrow = length( unique(.dat$facetVar1) ) ) 
  }
  
  
  # ~ Set Y-axis breaks ----------
  # other outcomes follow rules or can just use default axis breaks
  # y.breaks are only still null if none of the above applied
  if ( is.null(.y.breaks) ) {
    # set default breaks
    if ( grepl(pattern = "Cover", Yname) ){
      y.breaks = seq(0, 1, .1)
      
    } else {
      # otherwise keep the default limits from ggplot
      y.breaks = ggplot_build(p)$layout$panel_params[[1]]$y$breaks
    }
  } # end "if ( is.null(.y.breaks) )"
  
  
  # if user provided their own y.breaks
  if ( !is.null(.y.breaks) ) {
    y.breaks = .y.breaks
  }
  
  
  # use coord_cartesian so that lines/points that go outside limits look cut off
  #  rather than completely omitted
  p = p + coord_cartesian( ylim = c( min(y.breaks), max(y.breaks) ) ) +
    scale_y_continuous( breaks = y.breaks )
  
  if ( .writePlot == TRUE ) {
    my_ggsave( name = paste(Yname, "_plot.pdf", sep=""),
               .plot = p,
               .width = 10,
               .height = 10,
               .results.dir = .results.dir,
               .overleaf.dir = overleaf.dir.figs )
  }
  
  return(p)
  
}


# major fn for making the main-text figs (1 row per outcome, 3 outcomes per figure)
#  AND supplement figs (1 entire plot per outcome)
# important: to show some outcomes ONLY in supp (e.g., MhatTestReject), you should have global vars like this:
# YnamesMain = c("MhatBias", "MhatCover", "MhatWidth")
# YnamesSupp = c("MhatBias", "MhatCover", "MhatWidth",
#                "MhatTestReject")
# the outcome(s) that only go in supp should be at the end of list
sim_plot_multiple_outcomes = function(.hack,
                                      .y.breaks = NULL,
                                      .ggtitle = "",
                                      .local.results.dir = NA) {
  
  .dat = agg
  # when using different levels of t2w:
  # .dat$facetVar = paste( "t2a=", .dat$t2a, "; t2w=", .dat$t2w, sep = "")
  #.dat$facetVar = paste( "t2a = ", .dat$t2a, sep = "")
  .dat$facetVar = sqrt(.dat$t2a)
  
  .dat = .dat %>% filter(method.pretty %in% method.keepers &
                           Mu == 0.5 &
                           true.sei.expr == "0.02 + rexp(n = 1, rate = 3)" &
                           hack == .hack )
  
  
  # force ordering of methods to match 3_analyze_lodder.R
  correct.order = c("Uncorrected",
                    "SM",
                    "Unhacked only",
                    
                    "RTMA",
                    "MAN",
                    "SMKH")
  
  
  .dat$method.pretty = factor(.dat$method.pretty, levels = rev(correct.order))
  
  
  # ~~ Make plot for each outcome in YNamesMain ------------
  plotList = list()
  plotListSupp = list()
  
  for ( .Yname in YnamesSupp ) {
    
    i = which(YnamesSupp == .Yname)
    
    .dat$Y = .dat[[.Yname]]
    
    # ~~ Set ggplot color palette ----
    # to see all palettes:
    # par(mar=c(3,4,2,2))
    # display.brewer.all()
    
    # # SAVE - this is if you want to get colors automatically
    # n.colors.needed = length(unique(.dat$method.pretty))
    # .colors = brewer.pal(n = n.colors.needed, name = "Dark2")
    # if( length(.colors) > n.colors.needed ) .colors = .colors[1:n.colors.needed]
    # # this maps the colors onto levels of the factor
    # names(.colors) = levels( factor(.dat$method.pretty) )
    # # highlight certain methods
    # .colors[ names(.colors) == "RTMA" ] = "red"
    
    
    # set color palette to match 3_analyze_lodder.R
    .colors = c(#SMKH = "#1B9E77",
      MAN = "#ff9900",
      RTMA = "red",
      `Unhacked only` = "#3399ff",
      SM = "#00cc00",
      Uncorrected = "black")
    
    myColorScale = scale_colour_manual(values = .colors)
    
    # ~~ Set ggplot linetype scale ----
    # by default, dotted lines
    # but use solid lines for new proposed methods
    .lty = rep("dashed", nuni(.dat$method.pretty))
    names(.lty) = names(.colors)
    
    newMethods = c(#"SMKH",
      "MAN",
      "RTMA")
    
    .lty[ names(.lty) %in% newMethods ] = "solid"
    
    myLtyScale = scale_linetype_manual(values = .lty)
    
    # ~~ Set axis titles ---------
    
    # only label x-axis in last plot since they'll be combined
    if ( .Yname == YnamesMain[ length(YnamesMain) ] ) {
      .xlab = "Number of published nonaffirmative results"
    } else {
      .xlab = ""
    }
    
    
    .ylab = .Yname
    if ( .Yname == "MhatBias" ) .ylab = "Bias"
    if ( .Yname == "MhatCover" ) .ylab = "CI coverage"
    if ( .Yname == "MhatWidth" ) .ylab = "CI width"
    if ( .Yname == "MhatTestReject" ) .ylab = "Power"
    
    # ~ Make base plot ----------
    p = ggplot( data = .dat,
                aes( x = k.pub.nonaffirm,
                     y = Y,
                     color = method.pretty,
                     lty = method.pretty) ) +
      
      #geom_point() +
      
      geom_line() +
      
      # manually provided colors
      myColorScale +
      
      # manually provided linetypes
      myLtyScale +
      
      # base_size controls all text sizes; default is 11
      # https://ggplot2.tidyverse.org/reference/ggtheme.html
      theme_bw(base_size = 20) +
      
      #scale_x_log10( breaks = unique(.dp$n) )
      # use only some values
      #scale_x_log10( breaks = c(500, 1000) ) +
      
      xlab(.xlab) +
      scale_x_continuous( breaks = c(10, 40, 70, 100) ) +
      
      ylab(.ylab) +
      guides( color = guide_legend(title = "Method") ) +
      theme_bw() +
      theme( text = element_text(face = "bold"),
             
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()
             # reduce whitespace for combined plot
             # https://github.com/wilkelab/cowplot/issues/31
             #plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    
    # ~ Add reference lines ----------
    if ( str_contains(x = .Yname, pattern = "Cover") ) {
      p = p + geom_hline( yintercept = 0.95,
                          lty = 1,
                          color = "gray" ) 
      
    }
    
    if ( str_contains(x = .Yname, pattern = "Bias") ) {
      p = p + geom_hline( yintercept = 0,
                          lty = 1,
                          color = "gray" ) 
      
    }
    
    # ~ Add facetting ----------
    # this block needs to be after adding geom_hlines so that the lines obey the facetting
    p = p + facet_wrap( ~ facetVar,
                        ncol = length( unique(.dat$facetVar) ),
                        labeller = label_bquote( cols = tau[a] ~ "=" ~ .(facetVar) ) ) 
    
    
    # ~ Set Y-axis breaks ----------
    # other outcomes follow rules or can just use default axis breaks
    # y.breaks are only still null if none of the above applied
    if ( is.null(.y.breaks) ) {
      # set default breaks
      if ( grepl(pattern = "Cover", .Yname) ){
        y.breaks = seq(0, 1, .1)
        # possibly expanded breaks for supplement version
        y.breaks.supp = seq(0, 1, .1)
        
      } else if ( grepl(pattern = "Bias", .Yname) ){
        y.breaks = seq(-0.4, 0.5, .1)
        y.breaks.supp = seq(-1, 1, 0.1)
        
      } else if ( grepl(pattern = "Width", .Yname) ){
        y.breaks = seq(0, 2, .25)
        y.breaks.supp = seq(0, 3, 0.25)
        
      } else if ( grepl(pattern = "Reject", .Yname) ){
        y.breaks = seq(0, 1, .1)
        y.breaks.supp = y.breaks
        
      } else {
        # otherwise keep the default limits from ggplot
        y.breaks = ggplot_build(p)$layout$panel_params[[1]]$y$breaks
      }
    } # end "if ( is.null(.y.breaks) )"
    
    
    # if user provided their own y.breaks
    if ( !is.null(.y.breaks) ) {
      y.breaks = .y.breaks
    }
    
    
    # use coord_cartesian so that lines/points that go outside limits look cut off
    #  rather than completely omitted
    p = p + coord_cartesian( ylim = c( min(y.breaks), max(y.breaks) ) ) +
      scale_y_continuous( breaks = y.breaks )
    
    
    
    # ~ Handle ggtitle ----------------
    
    # only show title in the first plot since they'll be combined
    if ( i == 1 ) {
      p = p + ggtitle(.ggtitle)
    }
    
    # ~ Handle legend ----------------
    # combine legends into one
    p = p + labs(color  = "Method", linetype = "Method")
    
    # only show legend in the last plot since they'll be combined
    if ( .Yname == YnamesMain[ length(YnamesMain) ] ) {
      p = p + theme(legend.position = "bottom") +
        # fix order of legend items
        guides(colour = guide_legend(reverse = TRUE),
               linetype = guide_legend(reverse = TRUE) )
    } else {
      p = p + theme(legend.position = "none")
    }
    
    # some outcomes are being plotted only for supp
    #  don't include those in the plot list for the combined figure
    
    if ( .Yname %in% YnamesMain ) plotList[[i]] = p
    
    
    
    # ~ Modify plot for supplement --------------
    
    # expand axis limits and have 2 rows per outcome
    p = suppressMessages( p + coord_cartesian( ylim = c( min(y.breaks.supp), max(y.breaks.supp) ) ) +
                            scale_y_continuous( breaks = y.breaks.supp) +
                            
                            facet_wrap( ~ facetVar,
                                        ncol = 2,  # fix to 2 becuase will take up whole page
                                        labeller = label_bquote( cols = tau[a] ~ "=" ~ .(facetVar) ) ) + 
                            
                            # theme_bw(base_size=21) +
                            # theme(text = element_text(face = "bold"),
                            #       legend.position = "bottom") +
                            
                            # always show legend
                            theme(legend.position = "bottom") +
                            
                            # always show title
                            ggtitle(.ggtitle) +
                            # fix order of legend items
                            guides(colour = guide_legend(reverse = TRUE),
                                   linetype = guide_legend(reverse = TRUE) ) ) 
    
    
    plotListSupp[[i]] = p 
    
  }  # end "for ( Y in YnamesMain )"
  
  # ~~ Nicely arrange plots as columns ------------
  
  # give extra space to last one to accommodate y-axis label
  nOutcomes = length(YnamesMain)
  # if nOutcomes = 4, use 1.5 in last slot here
  rel.heights = c(rep(1, nOutcomes-1), 1.3)
  pCombined = cowplot::plot_grid(plotlist = plotList,
                                 nrow = nOutcomes,
                                 rel_heights = rel.heights)
  
  
  # ~~ Write combined plot to Overleaf and local dir ------------
  name = paste( tolower(.hack),
                "_Mu0.5.pdf",
                sep = "" )
  
  if ( overwrite.res == TRUE ) {
    my_ggsave( name = name,
               .plot = pCombined,
               .width = 8,
               .height = 11,
               .results.dir = paste(.local.results.dir, "Simple plots in main text", sep = "/"),
               .overleaf.dir = overleaf.dir.figs )
  } else {
    message("\n\nNot writing the plot to local dir or Overleaf because overwrite.res = FALSE")
  }
  
  # ~~ Save each individual supplement plot ------------
  
  if ( overwrite.res == TRUE ) {
    
    for ( .Yname in YnamesSupp ) {
      name = paste( "supp_",
                    tolower(.hack),
                    "_",
                    tolower(.Yname),
                    "_Mu0.5.pdf",
                    sep = "" )
      my_ggsave( name = name,
                 .plot = plotListSupp[[ which(YnamesSupp == .Yname) ]],
                 .width = 8,
                 .height = 11,
                 .results.dir = paste(.local.results.dir, "Simple plots in Supplement", sep = "/"),
                 .overleaf.dir = overleaf.dir.figs )
      
    }
    
  } else {
    message("\n\nNot writing the plot to local dir or Overleaf because overwrite.res = FALSE")
  }
  
  # ~~ Print plot and return the list --------------
  pCombined
  return(plotListSupp)
} 


# SAPH-SPECIFIC SMALL HELPERS ----------------------------------

# initialize global variables that describe estimate and outcome names, etc.
init_var_names = function() {
  
  ### Names of statistical metrics ###
  # used later to create plots and tables, but needed to check var types 
  #  upon reading in data
  estNames <<- c("Mhat", "Shat")
  
  # blank entry is to get Mhat itself, which is useful for 
  #  looking at whether MAON>0
  mainYNames <<- c("Bias", "", "RMSE", "Cover", "Width", "EmpSE")
  
  otherYNames <<- c("EstFail", "CIFail", "RhatGt1.01", "RhatGt1.05")
  
  # these ones don't fit in nicely because the "Mhat" is in the middle of string
  #"OptimxPropAgreeConvergersMhatWinner", "OptimxNAgreeOfConvergersMhatWinner"
  MhatMainYNames <<- paste( "Mhat", c(mainYNames), sep = "" )
  MhatYNames <<- c( paste( "Mhat", c(mainYNames, otherYNames), sep = "" ) )
  #"OptimxPropAgreeConvergersMhatWinner", "OptimxNAgreeOfConvergersMhatWinner" )
  
  
  ### Names of parameter variables ###
  # figure out which scen params were actually manipulated
  # **assumes that k.pub.nonaffirm is always the first scen parameter variable
  ( param.vars <<- names(agg)[ which( names(agg) == "k.pub.nonaffirm" ) : which( names(agg) == "method" ) ] )
  
  
  # how many levels does each param var have in dataset?
  ( n.levels <<- agg %>% dplyr::select(param.vars) %>%
      summarise_all( function(x) nuni(x) ) )
  
  ( param.vars.manip <<- names(n.levels)[ n.levels > 1 ] )
  
  
  # eliminate redundant ones
  if ( "t2a" %in% param.vars.manip ) param.vars.manip <<- drop_vec_elements( param.vars.manip, c("S", "V") )
  
  
  ( param.vars.manip2 <<- drop_vec_elements(param.vars.manip, "method") )
  
  cat( paste("\n\nManipulated parameter vars: ",
             paste(param.vars.manip2, collapse= ", ") ) )
  
}

# GENERIC SMALL HELPERS -------------------------------------------------------------

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}

# quick median with NAs removed
medNA = function(x){
  median(x, na.rm = TRUE)
}


# quick median with NAs removed and 10th and 90th percentiles
medNA_pctiles = function(x){
  paste( round( median(x, na.rm = TRUE), 2 ),
         " (",
         round( quantile(x, probs = 0.10, na.rm = TRUE), 2 ),
         ", ",
         round( quantile(x, probs = 0.90, na.rm = TRUE), 2 ),
         ")",
         sep = "" )
}



names_with = function(.dat, .pattern) {
  names(.dat)[ grepl(pattern = .pattern, x = names(.dat) ) ]
}


# one or both dirs can be NA
my_ggsave = function(name,
                     .plot = last_plot(),
                     .width,
                     .height,
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir) {
  
  dirs = c(.results.dir, .overleaf.dir)
  dirIsNA = sapply(dirs, is.na)
  validDirs = dirs[ !dirIsNA ]
  
  
  for ( dir in validDirs ) {
    setwd(dir)
    ggsave( name,
            plot = .plot,
            width = .width,
            height = .height,
            device = "pdf" )
  }
}

# drop elements from vector by their values
drop_vec_elements = function(x, 
                             values.to.drop) {
  x[ !(x %in% values.to.drop)]
}


# sort.Yname: UNQUOTED name of performance variable to sort on
# keepers: vars to retain in the dataset
sort_agg = function( sort.Yname,
                     desc = TRUE,
                     keepers = c("scen.name", param.vars.manip, MhatMainYNames) ) {
  
  if ( desc == TRUE ) {
    agg %>% select(keepers) %>%
      arrange( desc( {{sort.Yname}} ) ) %>%
      mutate_if( is.numeric, function(x) round(x, 2) )
  } else {
    agg %>% select(keepers) %>%
      arrange( {{sort.Yname}} ) %>%
      mutate_if( is.numeric, function(x) round(x, 2) )
  }
  
}

# quickly look at results when running locally
srr = function() {
  
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




# read/write intermediate work
write_interm = function(x, filename){
  setwd(prepped.data.dir)
  #setwd("Intermediate work")
  write.csv(x, filename)
}

read_interm = function(filename){
  setwd(prepped.data.dir)
  #setwd("Intermediate work")
  read.csv(filename)
}

# like View(), but opens the extra tab if global var useView = TRUE
View2 = function(x){
  if ( useView == TRUE ) View(x) 
}

# quick length(unique) equivalent
uni = function(x){
  length(unique(x))
}


# return strings containing anything in pattern vector
stringsWith = function(pattern, x){
  # make regex expression 
  patterns = paste(pattern, collapse="|")
  x[ grepl(pattern = patterns, x = x)]
}
# stringsWith( pattern = c("dog", "cat"),
#  x = c("dogcat", "horse", "cat", "lion") )


# return indices of strings containing anything in pattern vector
whichStrings = function(pattern, x){
  patterns = paste(pattern, collapse="|")
  grepl(pattern = pattern, x = x)
}

# stands for "wipe results"
wr = function(){
  #setwd(results.dir)
  #if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  setwd(overleaf.dir.nums)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
}

# stands for "view results"
vr = function(){
  setwd(overleaf.dir.nums)
  View( read.csv("stats_for_paper.csv") )
}


# make a string for estimate and CI
stat_CI = function(est, lo, hi){
  paste( est, " [", lo, ", ", hi, "]", sep = "" )
}
# stat_CI( c(.5, -.1), c(.3, -.2), c(.7, .0) )


# return percent true for 0/1 variable, counting NA as own category
percTRUE_incl_NA = function(x) {
  prop.table( table(x, useNA = "ifany") )[2]
}

# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
# expects "study" to be a global var
# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
# expects "study" to be a global var
update_result_csv = function( name,
                              .section = NA,
                              .results.dir = NULL,
                              .overleaf.dir = NULL,
                              value = NA,
                              print = FALSE ) {
  
  # if either is NULL, it just won't be included in this vector
  dirs = c(.results.dir, .overleaf.dir)
  
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(.section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  for (.dir in dirs) {
    
    setwd(.dir)
    
    if ( "stats_for_paper.csv" %in% list.files() ) {
      res <<- read.csv( "stats_for_paper.csv",
                        stringsAsFactors = FALSE,
                        colClasses = rep("character", 3 ) )
      
      # if this entry is already in the results file, overwrite the
      #  old one
      if ( all(name %in% res$name) ) res[ res$name %in% name, ] <<- new.rows
      else res <<- rbind(res, new.rows)
    }
    
    if ( ! "stats_for_paper.csv" %in% list.files() ) {
      res <<- new.rows
    }
    
    write.csv( res, 
               "stats_for_paper.csv",
               row.names = FALSE,
               quote = FALSE )
    
  }  # end "for (.dir in dirs)"
  
  
  if ( print == TRUE ) {
    View(res.overleaf)
  }
  
}


quick_ci = function( est, var ) {
  c( est - qnorm(.975) * sqrt(var),
     est + qnorm(.975) * sqrt(var) )
}

quick_pval = function( est, var ) {
  2 * ( 1 - pnorm( abs( est / sqrt(var) ) ) )
}


# # STRAIGHT FROM MRM:
# # summarize performance metrics given a dataset (dat) that is already scenario-aggregated
# #  looks for all variables with "Bias" or "Cover" in their names and takes their means
# # description: description of the row to make a nice table
# # selectVars: "Phat", "Diff" (by default looks for a global variable by this name)
# my_summarise = function(dat,
#                         description = NA,
#                         .selectVars = selectVars,
#                         badCoverageCutoff = 0.85,
#                         badWidthCutoff = 0.90,
#                         averagefn = "mean"
# ){
#   
#   # variables whose average should be taken
#   meanVars = c( namesWith(pattern = "Bias", dat = dat), 
#                 namesWith(pattern = "AbsErr", dat = dat),
#                 namesWith(pattern = "Cover", dat = dat),
#                 namesWith(pattern = "Width", dat = dat) )
#   
#   if (.selectVars == "Phat") meanVars = meanVars[ !grepl(pattern = "Diff", x = meanVars) ]
#   
#   if (.selectVars == "Diff") meanVars = meanVars[ !grepl(pattern = "Phat", x = meanVars) ]
#   
#   
#   if (averagefn == "mean") avgfun = function(x) meanNA(x)
#   if (averagefn == "median") avgfun = function(x) medNA(x)
#   if (averagefn == "median.pctiles") avgfun = function(x) medNA_pctiles(x)
#   
#   # make a one-row summary
#   tab = dat %>% 
#     summarise_at( .vars = meanVars, 
#                   function(x) avgfun(x) )
#   
#   # proportion of SCENARIOS with bad MEAN coverage
#   tab = tab %>% add_column(BadPhatCover = mean(dat$CoverPhat < badCoverageCutoff),
#                            BadDiffCover = mean(dat$CoverDiff < badCoverageCutoff))
#   if (.selectVars == "Phat") tab = tab %>% select(-BadDiffCover)
#   if (.selectVars == "Diff") tab = tab %>% select(-BadPhatCover)
#   
#   # proportion of SCENARIOS with bad median CI width
#   tab = tab %>% add_column(BadPhatWidth = mean(dat$PhatCIWidth > badWidthCutoff),
#                            BadDiffWidth = mean(dat$DiffCIWidth > badWidthCutoff) )
#   if (.selectVars == "Phat") tab = tab %>% select(-BadDiffWidth)
#   if (.selectVars == "Diff") tab = tab %>% select(-BadPhatWidth)
#   
#   # only round selected columns if we have strings
#   #  (e.g., median with percentiles)
#   if ( averagefn == "median.pctiles" ){
#     
#     if (.selectVars == "Phat"){
#       tab$BadPhatCover = round(tab$BadPhatCover, 2)
#       tab$BadPhatWidth = round(tab$BadPhatWidth, 2)
#     }
#     
#     if (.selectVars == "Diff"){
#       tab$BadDiffCover = round(tab$BadDiffCover, 2)
#       tab$BadDiffWidth = round(tab$BadDiffWidth, 2)
#     } 
#     
#   } else {  # otherwise round all columns
#     tab = round( tab, 2 )
#   }
#   
#   tab = tab %>% add_column(n.scens = nrow(dat), .before = 1 )
#   
#   if ( !is.na(description) ) tab = tab %>% add_column(Scenarios = description, .before = 1)
#   return(tab)
# }



