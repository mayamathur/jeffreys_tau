

# NOTES ----------------------------------------------------

# Only results from sim.env = stefan are in paper per reviewers' comments.
#  Analyses with sim.env = mathur retained here for completeness.


# PRELIMINARIES ----------------------------------------------------

#rm(list=ls())

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
# renv::restore()

# data-wrangling packages
library(here)
library(plotly)  # must be BEFORE dplyr or else plotly::select will take over
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(tidyverse)
library(fastDummies)
library(xlsx)
# meta-analysis packages
library(metafor)
library(robumeta)
# other
library(xtable)
library(testthat)
library(Deriv)
library(mosaic)
library(hpa)
library(pracma)
library(truncnorm)
library(tmvtnorm)
library(RColorBrewer)
library(sjmisc)

# prevent masking
select = dplyr::select

# run this only if you want to update the R environment specs
# setwd(here())
# renv::snapshot()

# ~~ User-specified global vars -------------------------
# no sci notation
options(scipen=999)

stitch.from.scratch = TRUE

# control which results should be redone and/or overwritten
# but note that not all fns respect this setting
overwrite.res = TRUE


# ~~ Set directories -------------------------
code.dir = here()

( data.dir = str_replace( string = here(),
                          pattern = "Code \\(git\\)",
                          replacement = "Results/Working dataset") )

( results.dir = str_replace( string = here(),
                             pattern = "Code \\(git\\)",
                             replacement = "Results/Working results") )

# check that they're specified correctly
setwd(data.dir)
setwd(results.dir)

# below is the only absolute path
# write results directly to directory containing TeX manuscript in Overleaf so stats can be piped directly into text
# this is an absolute path because it must live in Dropbox, outside the project directory, in order to sync with Overleaf
# to reproduce results, just set this to any directory on your local machine
# results will be written to a csv file in that location
overleaf.dir.figs = results.dir


setwd(code.dir)
source("analyze_sims_helper_JTE.R")


Ynames = c("ShatAbsBias", "ShatCover", "ShatRMSE", "ShatWidth",
           "MhatAbsBias", "MhatCover", "MhatRMSE", "MhatWidth")

param.vars = c("true.dist", "true.sei.expr", "k.pub", "Mu", "t2a")




# ~~ Get agg data -------------------------

# if only analyzing a single sim environment (no merging):
setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agg)


# drop any "NA" methods (i.e., ones that didn't get labeled in wrangle_agg_local)
agg = agg %>% filter( method.pretty != "" )
table(agg$method.pretty)

# look at number of actual sim reps
table(agg$sim.reps.actual)



# ~~ List variable names -------------------------

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()


# ~~ Make data subsets -------------------------

# realistically small metas only
aggs = agg %>% filter( k.pub <= 20 )



# ~~ Convergence stats by method -------------------------

summary(aggs$MhatEstConverge)
summary(aggs$MhatCIFail)

aggs %>% group_by(method)


# convergence rates
t = aggs %>% group_by(method) %>%
  summarise( mean(1-MhatEstFail), 
             min(1-MhatEstFail),
             
             mean(1-MhatCIFail),
             min(1-MhatCIFail) )

View(t)

# PLOTS -------------------------------------------------

.true.dist = "norm"
.true.sei.expr = "0.02 + rexp(n = 1, rate = 3)"
.Mu = 0.5
.estimate = "Shat"
.local.results.dir = results.dir




sim_plot_multiple_outcomes(.agg = agg,
                           
                           .true.dist = .true.dist,
                           .true.sei.expr = .true.sei.expr,
                           .Mu = .Mu,
                           .estimate = .estimate, 
                           
                           .y.breaks = NULL,
                           .ggtitle = "",
                           .local.results.dir = .local.results.dir)



# QUICK AND SIMPLE SUBSET ANALYSIS -------------------------------------------------


t = aggo %>%
  filter(k.pub <= 20) %>%
  #filter(true.dist == "expo") %>%
  #filter(true.sei.expr == "0.02 + rexp(n = 1, rate = 1)") %>%
  filter(true.sei.expr != "0.3") %>%
  group_by(method) %>%
  summarise( meanNA(ShatBias),
             meanNA(ShatCover),
             meanNA(ShatRMSE),
             meanNA(ShatWidth),
             
             meanNA(MhatBias),
             meanNA(MhatCover),
             meanNA(MhatRMSE),
             meanNA(MhatWidth) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))

View(t)

# ******** WINNER TABLES -------------------------

make_both_winner_tables(.agg = agg)

# small metas
make_both_winner_tables(.agg = agg %>% filter(k.pub <= 20))
make_both_winner_tables(.agg = agg %>% filter(k.pub == 5))

# t2a
make_both_winner_tables(.agg = agg %>% filter(t2a == 0.0025 & k.pub < 20))
make_both_winner_tables(.agg = agg %>% filter(t2a == 1 & k.pub < 20))



# PLOTS -------------------------------------------------


# EXPLORE PREDICTORS OF PERFORMANCE  -------------------------------------------------

# not sure how useful this is



# possible figures:
# x = k.pub
# colors = methods
# row panels = Yname
# column panels = t2a

# average over: true.dist, true.sei.expr, Mu?

regressions.from.scratch = TRUE




performance_regressions = function(.agg,
                                   Ynames = Ynames,
                                   covariates = param.vars ) {
  
  
  for (i in Ynames) {
    
    # # TEST
    # .agg = agg %>% filter(method == "jeffreys-pmed")
    # Ynames = Ynames
    # covariates = param.vars
    # i = "ShatCover"
    
    string = paste( i, "~", paste(covariates, collapse = "+"), sep="" )
    mod = lm( eval(parse(text=string)),
              data = .agg )
    
    coefs = coef(mod)[-1]
    
    pvals = summary(mod)$coefficients[,"Pr(>|t|)"]
    pvals = pvals[-1]  # remove intercept
    
    # which vars are good vs. bad for the outcome?
    # flip coeff signs so that positive means it improves the outcome
    if ( grepl(pattern = "AbsBias", x = i) | grepl(pattern = "Width", x = i) | grepl(pattern = "RMSE", x = i) ) coefs = -coefs
    good = names( coefs[ coefs > 0 & pvals < 0.01 ] )
    bad = names( coefs[ coefs < 0 & pvals < 0.01 ] )
    
    good = paste(good, collapse = ", ")
    bad = paste(bad, collapse = ", ")
    
    newRow = data.frame( outcome = i,
                         good = good, 
                         bad = bad )
    
    if (i==Ynames[1]) res = newRow else res = rbind(res, newRow)
    
    cat( paste("\n\n*******************", toupper(i), " PREDICTORS*******************\n" ) )
    print(summary(mod))
    
  }  # end loop over outcomes
  
  
  # look at results
  res
  
  # clean up string formatting
  # res2 = res %>% mutate_at( vars(good, bad), my_recode )
  
}


( t1 = performance_regressions(.agg = agg %>% filter(method == "jeffreys-pmed"),
                               Ynames = Ynames,
                               covariates = param.vars) )


( t2 = performance_regressions(.agg = agg %>% filter(method == "EB"),
                               Ynames = Ynames,
                               covariates = param.vars) )

