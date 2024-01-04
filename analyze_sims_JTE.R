

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

# dataset prepped by prep_sims_SAPH.R
data.dir = str_replace( string = here(),
                           pattern = "Code \\(git\\)/Sherlock code",
                           replacement = "Simulation results/*2023-06-21 aggregated simulations (as in published paper)" )

# check that they're specified correctly
setwd(data.dir)
setwd(results.dir)

# below is the only absolute path
# write results directly to directory containing TeX manuscript in Overleaf so stats can be piped directly into text
# this is an absolute path because it must live in Dropbox, outside the project directory, in order to sync with Overleaf
# to reproduce results, just set this to any directory on your local machine
# results will be written to a csv file in that location
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/P-hacking (SAPH)/figures_SAPH/sims"


setwd(code.dir)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")




# ~~ Get agg data -------------------------

# if only analyzing a single sim environment (no merging):
setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agg)


# drop any "NA" methods (i.e., ones that didn't get labeled in wrangle_agg_local)
agg = agg %>% filter( method.pretty != "" )

# look at number of actual sim reps
table(agg$sim.reps.actual)


#@important: stefan and mathur share scen.names, so relabel them
agg$scen.name2 = paste(agg$sim.env, agg$scen.name)
# stefan: 40 scens, mathur: 48
expect_equal( nuni(agg$scen.name2), 40 + 48)
 


# ~~ List variable names -------------------------

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()


# ~~ Make data subsets -------------------------

# stefan only
aggs = agg %>% filter(sim.env == "stefan")
# mathur only
aggm = agg %>% filter(sim.env == "mathur")



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



# ******** WINNER TABLES -------------------------

# ~ Winner tables in paper ------------------------------

# can toggle output of fn below by changing the default arg of 
#  make_winner_table between display = "dataframe" (easy viewing)
#  and display = "xtable" (Overleaf)



# ~~~ Stefan ------------------------------
# all scenarios
make_both_winner_tables(.agg = aggs)

# by k
make_both_winner_tables(.agg = aggs %>% filter(k.pub.nonaffirm == 10) )


# by two- vs. one-tailed selection
make_both_winner_tables(.agg = aggs %>% filter(alternative.stefan == "two.sided") )
make_both_winner_tables(.agg = aggs %>% filter(alternative.stefan == "greater") )

# by strategy (favor-first vs. favor-best)
make_both_winner_tables(.agg = aggs %>% filter(strategy.stefan == "firstsig") )
make_both_winner_tables(.agg = aggs %>% filter(strategy.stefan == "smallest") )



# ~ Winner tables not in paper per reviewers: Mathur sim environment ------------------------------

# can toggle output of fn below by changing the default arg of 
#  make_winner_table between display = "dataframe" (easy viewing)
#  and display = "xtable" (Overleaf)

# ~~~ Mathur  ------------------------------

# 1: scenarios where stringent overall selection holds (not "favor-best-affirm-wch", "affirm")
make_both_winner_tables(.agg = aggm %>% filter(rtma.misspec == FALSE))

# 2: all scenarios
make_both_winner_tables(.agg = aggm)
# scenarios with "favor-best-affirm-wch" or "affirm"
make_both_winner_tables(.agg = aggm %>% filter(rtma.misspec == TRUE))



