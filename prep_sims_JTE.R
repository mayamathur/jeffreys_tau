
# NOTES ----------------------------------------------------

# See notes in doParallel_JTE.R about the two simulation batches. 
# Make sure you set the global variable sim_set below depending on the batch to be analyzed.


# PRELIMINARIES ----------------------------------------------------

#  rm(list=ls())

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
library(tableone)

# prevent masking
select = dplyr::select

# run this only if you want to update the R environment specs
# renv::snapshot()

# no sci notation
options(scipen=999)

# ~~ User-specified global vars -------------------------

# are we running the main analysis, or the supplementary bootstrap analysis?
# this avoids results in stats_for_paper.csv
sim_set = "boot"
#sim_set = "main"
message("\n\n***** Setting sim_set = ", sim_set)


# ~~ Set directories -------------------------
code.dir = here()

# "official" directory names:
if ( sim_set == "main" ) {
  ( data.dir = str_replace( string = here(),
                            pattern = "Code \\(git\\)",
                            replacement = "Results/*2024-02-26 - collect pmed, mean from bayesmeta (as in RSM_0)/Datasets") )
  
  ( results.dir = str_replace( string = here(),
                               pattern = "Code \\(git\\)",
                               replacement = "Results/*2024-02-26 - collect pmed, mean from bayesmeta (as in RSM_0)/Results") )
}

if ( sim_set == "boot" ) {
  ( data.dir = str_replace( string = here(),
                            pattern = "Code \\(git\\)",
                            replacement = "Results/*2024-07-13 - add two types of boot (k=10 scens only)/Datasets") )
  
  ( results.dir = str_replace( string = here(),
                               pattern = "Code \\(git\\)",
                               replacement = "Results/*2024-07-13 - add two types of boot (k=10 scens only)/Results") )
}

# # generic directories (SAVE):
# ( data.dir = str_replace( string = here(),
#                           pattern = "Code \\(git\\)",
#                           replacement = "Results/Working dataset") )
# 
# ( results.dir = str_replace( string = here(),
#                              pattern = "Code \\(git\\)",
#                              replacement = "Results/Working results") )

# check that they're specified correctly
setwd( data.dir)
setwd(results.dir)


setwd(code.dir)
source("helper_JTE.R")
source("analyze_sims_helper_JTE.R")




# READ IN AGGREGATED DATA FROM CLUSTER -------------------------------------------------

# ~ Basic prep -------------------------------------------------
setwd(data.dir)
aggo = fread("aggo.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("aggo.csv")$mtime
nrow(aggo) / nuni(aggo$method)  # number of scens that are done; 3120 if sims are done
nuni(aggo$scen.name)

# add fancy variables for plotting, etc.
agg = wrangle_agg_local(aggo)
table(agg$method.pretty)

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names(.agg = agg)

agg = agg %>% filter(!is.na(scen.name))


# ~ Individual studies should be unbiased -------------------------------------------------

# look for scens where even the individual studies are biased for Mu
#  e.g., because of very rare binary Y with small N
summary( abs(agg$sancheck_mean_yi - agg$Mu) )

# flag scens where yi has bias > 0.05
agg$exclude_scen_biased_yi = abs(agg$sancheck_mean_yi - agg$Mu) > 0.05
mean(agg$exclude_scen_biased_yi)  # percent of scens

message( paste( "\n\n", round( 100 * meanNA(agg$exclude_scen_biased_yi) ), "% of scens had biased yi and will be removed", sep = " ") )


# expect no NAs here
mean(is.na(agg$sancheck_mean_yi))

# summarize scen params for these ones
# not surprisingly, the bad scens are exclusively binary Y, and mostly ones with N=40 or N ~ Unif(40,400)
#  though spread across a variety of p0 values
agg_bad = agg %>% filter(exclude_scen_biased_yi == TRUE)
CreateTableOne( dat = agg_bad,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2 )

table(agg_bad$N.expr )

# write dataset before excluding scens
setwd(data.dir)
fwrite(agg, "agg_including_scens_biased_yi.csv")


#***exclude these scens going forward
agg = agg %>% filter(exclude_scen_biased_yi == FALSE)
CreateTableOne( dat = agg,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2,
                strata = "Ytype")

nuni(agg$scen.name)


# Write prepped datasets -------------------------------------------------

setwd(data.dir)
fwrite(agg, "agg.csv")
fwrite(agg_bad, "agg_just_the_excluded_scens_biased_yi.csv")


# PREP ITERATE-LEVEL DATA FOR SCEN 1384  -------------------------------------------------

if ( sim_set == "main" ) {
  # this will be a little slow (1-2 min)
  setwd(data.dir)
  s2 = fread("long_results_job_1384.csv")
  
  expect_equal( 500, nrow(s2) / nuni(s2$method) )
  
  # make analysis vars
  s2 = s2 %>% rowwise() %>%
    mutate( CI_asy = (MHi - Mhat) / (Mhat - MLo),
            MhatBias = Mhat - Mu,
            MhatWidth = MHi - MLo,
            MhatCover = (MHi >= Mu & MLo <= Mu) )
  
  # recode variables
  s2$method.pretty.mu.inf = s2$method 
  s2$method.pretty.mu.inf[ s2$method == "bayesmeta-joint-shortest-margpmode" ] = "Jeffreys2-shortest" 
  s2$method.pretty.mu.inf[ s2$method == "bayesmeta-tau-shortest-margpmode" ] = "Jeffreys1-shortest"
  
  s2$method.pretty.mu.inf[ s2$method == "ML" ] = "ML-HKSJ"
  s2$method.pretty.mu.inf[ s2$method == "PM" ] = "PM-HKSJ"
  s2$method.pretty.mu.inf[ s2$method == "DL" ] = "DL-HKSJ"
  s2$method.pretty.mu.inf[ s2$method == "DL2" ] = "DL2-HKSJ"
  s2$method.pretty.mu.inf[ s2$method == "REML" ] = "REML-HKSJ"
  s2$method.pretty.mu.inf[ s2$method == "exact" ] = "Exact"
  
  
  fwrite( s2, "pretty_long_results_job_1384.csv" )
}

