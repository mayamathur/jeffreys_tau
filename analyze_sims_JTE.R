

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
library(tableone)

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


# ~~ Get agg data -------------------------

# if only analyzing a single sim environment (no merging):
setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agg) / nuni(agg$method)


# drop any "NA" methods (i.e., ones that didn't get labeled in wrangle_agg_local)
agg = agg %>% filter( method.pretty != "" )
table(agg$method.pretty)

# look at number of actual sim reps
table(agg$sim.reps.actual)



# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()


# summarize scen params
CreateTableOne( dat = agg,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2,
                strata = "Ytype" )


# ~~ Check runtimes of sbatch files -------------------------

# mean runtimes within scenarios - HOURS
summary(agg$doParallelSeconds/60^2) 

# 95th quantile of runtime within scens - HOURS
summary(agg$doParallelSecondsQ95/60^2) 



# ~~ Make data subsets -------------------------

# agg_save = agg  # in case you want to subset
# 
# # realistically small metas only
# aggs = agg %>% filter( k.pub <= 20 )



# ~~ Convergence stats by method -------------------------

summary(agg$MhatEstFail)
summary(agg$MhatCIFail)

summary(agg$ShatEstFail)
summary( agg$ShatCIFail[ agg$method != "robu" ] ) # robu doesn't even try to provide inference

# convergence rates
t = agg %>% group_by(method) %>%
  summarise( mean(1-MhatEstFail), 
             min(1-MhatEstFail),
             
             mean(1-MhatCIFail),
             min(1-MhatCIFail),
             
             mean(1-ShatEstFail), 
             min(1-ShatEstFail),
             
             mean(1-ShatCIFail),
             min(1-ShatCIFail) )

View(t)



# SANITY CHECKS ON DATA GENERATION -------------------------


# ~ Individual studies should be unbiased  -------------------------------------------------

# look for scens where even the individual studies are biased for Mu
#  e.g., because of very rare binary Y with small N
summary( abs(agg$sancheck_mean_yi - agg$Mu) )

ind = which( abs(agg$sancheck_mean_yi - agg$Mu) > 0.1 )
length(ind)/nrow(agg)  # percent of scens

# summarize scen params for these ones
CreateTableOne( dat = agg[ind,],
                vars = param.vars.manip2,
                factorVars = param.vars.manip2)

table(agg[ind,"N.pretty"] )
# not surprisingly, the bad scens are exclusively binary Y, and almost exclusively ones with N=40
#  though spread across a variety of p0 values


# now with a more stringent threshold for bias in yi
ind = which( abs(agg$sancheck_mean_yi - agg$Mu) > 0.05 )
ind
length(ind)/nrow(agg)  # percent of scens

# summarize scen params for these ones
CreateTableOne( dat = agg[ind,],
                vars = param.vars.manip2,
                factorVars = param.vars.manip2 )

table(agg[ind,"N.pretty"] )
# not surprisingly, the bad scens are exclusively binary Y, and almost exclusively ones with N=40
#  though spread across a variety of p0 values


#***exclude these scens going forward
agg = agg[-ind,]
CreateTableOne( dat = agg,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2,
                strata = "Ytype")

# ~ Other sanity checks -------------------------------------------------


namesWith(pattern = "sancheck_", agg)

summary( abs( agg$sancheck_mean_pY0 - agg$p0 ) )
summary( abs( agg$sancheck_mean_nY0 - agg$sancheck_mean_nY0_theory) )
summary( abs( agg$sancheck_mean_nY1 - agg$sancheck_mean_nY1_theory) )


# PLOTS -------------------------------------------------


#bm: edit this fn to aggregate over vars not used in facetting

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





# AUTO-FIND INTERESTING SCENS  -------------------------------------------------

# scens with meaningful differences across methods in AbsBias, RMSE, or Cover for Mhat or Shat
# among those, label the winning method? or assign points for winning each of these?

# basically no meaningful differences in: 

# t = agg %>% group_by(scen.name) %>%
#   summarise( MhatAbsBiasSD = sd(MhatAbsBias, na.rm = TRUE),
#              MhatAbsBiasSD = sd(MhatAbsBias, na.rm = TRUE))

# alternate: range instead of SD
t = agg %>% group_by(scen.name) %>%
  summarise( MhatAbsBiasRange = diff( range(MhatAbsBias, na.rm = TRUE) ),
             MhatCoverRange = diff( range(MhatCover, na.rm = TRUE) ),
             MhatRMSERange = diff( range(MhatRMSE, na.rm = TRUE) ),
             
             ShatAbsBiasRange = diff( range(ShatAbsBias, na.rm = TRUE) ),
             ShatCoverRange = diff( range(ShatCover, na.rm = TRUE) ),
             ShatRMSERange = diff( range(ShatRMSE, na.rm = TRUE) ) )

summary(t$MhatAbsBiasRange)
summary(t$MhatCoverRange)
summary(t$MhatRMSERange)

summary(t$ShatAbsBiasRange)
summary(t$ShatCoverRange)
summary(t$ShatRMSERange)

# *much more variability on Shat performance metrics than on Mu

# mark scens as important if methods vary on any of these characteristics
t$scen_important_Mhat = t$MhatCoverRange > 0.1 | t$MhatRMSERange > 0.25
t$scen_important_Shat = t$ShatAbsBiasRange > 0.1 | t$ShatCoverRange > 0.1 | t$ShatRMSERange > 0.25
t$scen_important = t$scen_important_Mhat | t$scen_important_Shat

mean(t$scen_important_Mhat)
mean(t$scen_important_Shat)
mean(t$scen_important_Mhat | t$scen_important_Shat)  # the Mhat important ones are mostly a subset of the Shat important ones

important_scens = t$scen.name[ t$scen_important_Mhat | t$scen_important_Shat ]

# **add the importance vars to agg
agg = left_join(x = agg,
                y = t %>% select(scen.name, scen_important),
                by = "scen.name")


# experiment with how to auto-choose the winner
t2 = agg %>% filter(scen.name == 30) 

( x = make_winner_table_col(.agg = agg %>% filter(scen.name == 29),
                            yName = "MhatCover",
                            methods = unique(agg$method) ) )

( x = make_winner_table(.agg = agg %>% filter(scen.name == 29),
                        summarise.fun.name = "median") )



#bm: in t, have a variable for which subset of the 6 variables were important
# e.g., MhatAbsBias and ShatCover
# then 
# does this make sense??


# ******** WINNER TABLES -------------------------


.Ytype = "cont-SMD"
#.Ytype = "bin-OR"

# create the base dataset from which to filter all winner tables
agg2 = agg %>% filter( scen_important == TRUE & Ytype == .Ytype )
#agg2 = agg %>% filter( scen_important == TRUE & Ytype == .Ytype & t2a > 0.0001 )
#agg2 = agg %>% filter( Ytype == .Ytype & t2a > 0.0001 )


dim(agg2)
# summarize scen params
CreateTableOne( dat = agg2[ !duplicated(agg2$scen.name) ],
                vars = param.vars.manip2,
                factorVars = param.vars.manip2 )


# overall
make_both_winner_tables(.agg = agg2)
# what's going on with Mhat coverage here?
t = agg2 %>% group_by(method, Mu) %>%
  summarise( medianNA(MhatCover),
             medianNA(MhatWidth),
             medianNA(Mhat),
             meanNA(MLo), 
             meanNA(MHi) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))

View(t)

# power
# wow, huge difference in power
make_both_winner_tables(.agg = agg2 %>% filter(Mu > 0),
                        .yNames = "MhatTestReject" )
# false-positive rate
make_both_winner_tables(.agg = agg2 %>% filter(Mu == 0),
                        .yNames = "MhatTestReject" )


# small metas
make_both_winner_tables(.agg = agg2 %>% filter(k.pub <= 20))
make_both_winner_tables(.agg = agg2 %>% filter(k.pub == 2))
make_both_winner_tables(.agg = agg2 %>% filter(k.pub == 100))

# t2a: definitely matters
make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.0001))  # very bad for jeffreys (esp. Shat Cover)
make_both_winner_tables(.agg = agg2 %>% filter(t2a > 0.0001))  # good for jeffreys
make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.01))
make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.04)) # good for Jeffreys


# stratified by distribution
# makes little difference for all methods
make_both_winner_tables(.agg = agg2 %>% filter(true.dist == "norm"))
make_both_winner_tables(.agg = agg2 %>% filter(true.dist == "expo"))

# effect of having equal sample sizes vs. uniform
make_both_winner_tables(.agg = agg2 %>% filter(N.pretty == "N = 40"))
make_both_winner_tables(.agg = agg2 %>% filter(N.pretty == "N = 400"))
make_both_winner_tables(.agg = agg2 %>% filter( N.pretty == "N ~ U(40, 400)" ))
make_both_winner_tables(.agg = agg2 %>% filter( N.pretty == "N ~ U(2000, 3000)" ))




# 2024-01-13 - PULL OUT A PROBLEM SCENARIO -------------------------------------------------

x = agg2 %>% filter(t2a == 0.0001 & ShatCover == 0)
unique(x$scen.name)

x2 = agg %>% filter(scen.name==105); View(x2)

# scen params only
x3 = x2[1, 1:11] 
library(constructive)
construct( as.data.frame(x3) )



# SANITY CHECK ON WINNER TABLES: QUICK AND SIMPLE SUBSET ANALYSIS -------------------------------------------------


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



# EXPLORE PREDICTORS OF PERFORMANCE  -------------------------------------------------

# not sure how useful this is since what we really care about is relative performance



# possible figures:
# x = k.pub
# colors = methods
# row panels = Yname
# column panels = t2a

# average over: true.dist, true.sei.expr, Mu?

regressions.from.scratch = TRUE

# remove p0 because causes error "contrasts can be applied only to factors with 2 or more levels" in lm
#  because p0 doesn't vary conditional on Ytype = cont
param.vars.manip3 = drop_vec_elements(param.vars.manip2, "p0")


( t1 = performance_regressions(.agg = agg %>% filter(method == "jeffreys-max-lp-iterate"),
                               Ynames = Ynames,
                               covariates = param.vars.manip3 ) )


( t2 = performance_regressions(.agg = agg %>% filter(method == "EB"),
                               Ynames = Ynames,
                               covariates = param.vars.manip3) )

