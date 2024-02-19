

# NOTES ----------------------------------------------------


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
library(constructive)

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

# should all the "View()" things be called?
# they open new RStudio tabs, so can be useful to skip these sanity checks
use.View = TRUE

# ~~ Set directories -------------------------
code.dir = here()

( data.dir = str_replace( string = here(),
                          pattern = "Code \\(git\\)",
                          replacement = "Results/*2024-01-31 - as in RSM_0/Datasets") )

( results.dir = str_replace( string = here(),
                             pattern = "Code \\(git\\)",
                             replacement = "Results/*2024-01-31 - as in RSM_0/Results") )

# # generic directories (SAVE):
# ( data.dir = str_replace( string = here(),
#                           pattern = "Code \\(git\\)",
#                           replacement = "Results/Working dataset") )
# 
# ( results.dir = str_replace( string = here(),
#                              pattern = "Code \\(git\\)",
#                              replacement = "Results/Working results") )

# check that they're specified correctly
setwd(data.dir)
setwd(results.dir)

# below are the only absolute paths
# write results directly to directory containing TeX manuscript in Overleaf so stats can be piped directly into text
# this is an absolute path because it must live in Dropbox, outside the project directory, in order to sync with Overleaf
# to reproduce results, just set this to any directory on your local machine
# results will be written to a csv file in that location
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/JTE (Jeffreys tau estimation) Overleaf/R_objects/figures"
overleaf.dir.stats = "/Users/mmathur/Dropbox/Apps/Overleaf/JTE (Jeffreys tau estimation) Overleaf/R_objects"


setwd(code.dir)
source("analyze_sims_helper_JTE.R")
source("helper_JTE.R")  # for lprior(), etc.




# ~~ Read datasets -------------------------

# if only analyzing a single sim environment (no merging):
setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agg) / nuni(agg$method)

# reorder factors for ggplot
agg$Ytype.pretty = factor(agg$Ytype.pretty, levels = c("Continuous Y", "Binary Y") )
# reorder them
agg$true.dist.pretty = factor( agg$true.dist.pretty, levels = c("Normal effects", "Exponential effects"))



# main analysis dataset
agg2 = agg %>% filter(k.pub <= 20)


# # drop any "NA" methods (i.e., ones that didn't get labeled in wrangle_agg_local)
# agg = agg %>% filter( method.pretty != "" )
# table(agg$method.pretty)

# check that all sim reps completed
expect_equal( unique(agg$sim.reps.actual), 500 )

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names(.agg = agg)

# # summarize scen params
# CreateTableOne( dat = agg,
#                 vars = param.vars.manip2,
#                 factorVars = param.vars.manip2,
#                 strata = "Ytype" )



# ~~ Check runtimes of sbatch files -------------------------

# mean runtimes within scenarios - HOURS
summary(agg$doParallelSeconds/60^2) 

# 95th quantile of runtime within scens - HOURS
summary(agg$doParallelSecondsQ95/60^2) 



# SANITY CHECKS -------------------------

# ~ Compare methods that should be similar or identical  -------------------------------------------------

# ~~ MLE-profile vs. ML-Wald-Qprofile  -------------------------------------------------

# CIs are quite different for these two
# shouldn't CIs for mu be the same for these two?
t = agg %>% filter(method.pretty %in% c("MLE-Wald-Qprofile", "MLE-profile")) %>%
  group_by(scen.name) %>%
  summarise( sd(Mhat),
             sd(Shat),
             sd(MLo))

summary(t$`sd(Mhat)`)
summary(t$`sd(Shat)`)
summary(t$`sd(MLo)`)  # CIs are NOT always the same


# ~~ metaLik should be equivalent to MLE-profile  -------------------------------------------------
# are they always the same?
# very close, but not exact
t = agg %>% filter(method.pretty %in% c("metaLik", "MLE-profile")) %>%
  group_by(scen.name) %>%
  summarise( sd(Mhat),
             sd(Shat),
             sd(MLo))

summary(t$`sd(Mhat)`)
summary(t$`sd(Shat)`)
summary(t$`sd(MLo)`)  

# their average performances are the same, though:
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" & method.pretty %in% c("metaLik", "MLE-profile") ) )
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" & method.pretty %in% c("metaLik", "MLE-profile") ) )


# ~~ Bayesian methods -------------------------------------------------

t = agg %>% filter(method %in% c("bayesmeta-joint-central", "bayesmeta-joint-shortest")) %>%
  group_by(scen.name) %>%
  summarise( sd(Mhat),
             sd(Shat),
             sd(MLo))

summary(t$`sd(Mhat)`)
summary(t$`sd(Shat)`)
summary(t$`sd(MLo)`) 

t = agg %>% filter(method %in% c("bayesmeta-tau-central", "bayesmeta-tau-shortest")) %>%
  group_by(scen.name) %>%
  summarise( sd(Mhat),
             sd(Shat),
             sd(MLo))

summary(t$`sd(Mhat)`)
summary(t$`sd(Shat)`)
summary(t$`sd(MLo)`) 


# ~ One-off stats for paper -------------------------------------------------

### Stats about scen parameters
# one row per scen only
first = agg[ !duplicated(agg$scen.name), ]

update_result_csv( name = "Num scens Ytype bin",
                   value = sum(first$Ytype == "bin-OR"),
                   print = TRUE )

update_result_csv( name = "Num scens Ytype cont",
                   value = sum(first$Ytype == "cont-SMD"),
                   print = TRUE )



### CI width comparisons
t = agg %>% filter(Ytype == "cont-SMD") %>%
  group_by(scen.name) %>%
  mutate( CI_ratio = min( MhatWidth[ method.pretty != "Jeffreys" ] ) / MhatWidth[ method.pretty == "Jeffreys" ] ) %>%
  filter( !duplicated(scen.name) )

expect_equal( nrow(t), nuni(agg$scen.name[agg$Ytype == "cont-SMD"]) )

#@definitely check these and the underlying CI_ratio calculation
update_result_csv( name = "Sims - Mean perc narrower Jeffreys vs winning other method - Ycont",
                   value = round( 100 * ( mean(t$CI_ratio) - 1 ) ),
                   print = TRUE )



# ~ Other sanity checks -------------------------------------------------

# sanity checks on data generation
namesWith(pattern = "sancheck_", agg)

summary( abs( agg$sancheck_mean_pY0 - agg$p0 ) )
summary( abs( agg$sancheck_mean_nY0 - agg$sancheck_mean_nY0_theory) )
summary( abs( agg$sancheck_mean_nY1 - agg$sancheck_mean_nY1_theory) )


# ~ Convergence stats by method -------------------------

# convergence rates
t = agg %>% group_by(method.pretty) %>%
  summarise( mean(1-MhatEstFail), 
             min(1-MhatEstFail),
             
             mean(1-MhatCIFail),
             min(1-MhatCIFail),
             
             mean(1-ShatEstFail), 
             min(1-ShatEstFail),
             
             mean(1-ShatCIFail),
             min(1-ShatCIFail) )

if (use.View == TRUE) View(t)


# Bayesian convergence metrics
( t = agg %>% filter(method.pretty == "Jeffreys") %>%
    summarise( mean(MhatRhatGt1.05),
               mean(ShatRhatGt1.05),
               max(MhatRhatGt1.05),
               max(ShatRhatGt1.05) ) )


update_result_csv( name = paste("Perc Jeffreys", names(t)),
                   value = round( 100 * t[1,], 2 ),
                   print = TRUE )


# WINNER TABLES -------------------------


dput(unique(agg$method))


# fewer BY methods
# methods_for_table = c(
#   "ML", "REML", "DL", "PM", "DL2", "exact","MLE-profile",
#   "bayesmeta-tau-central", 
#   "bayesmeta-tau-shortest",
#   "bayesmeta-joint-shortest",
#   "bayesmeta-joint-central")

# # create the base dataset from which to filter all winner tables
agg2 = agg %>% filter( k.pub <= 20 )


dim(agg2); nuni(agg2$scen.name)
# summarize scen params
CreateTableOne( dat = agg2[ !duplicated(agg2$scen.name) ],
                vars = param.vars.manip2,
                factorVars = param.vars.manip2 )

#  ~ Sanity checks -------------------------------------------------

# only look at methods that should be exactly the same

# very similar, but not exactly the same
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" &
                                                  method %in% c("bayesmeta-joint-shortest", "jeffreys-hdi") ) )

# **big difference here: bayesmeta is way shorter
# could it be that the posterior is multimodal?
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" &
                                                  method %in% c("bayesmeta-tau-shortest", "jeffreys-tau-hdi") ) )

# again extremely different! 
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" &
                                                  method %in% c("bayesmeta-tau-central", "jeffreys-tau-pmode") ) )


# ~ Overall  -------------------------------------------------

# if the tables don't have all the methods you want, adjust args in make_winner_table_col

# **MAIN TEXT TABLES 2-5
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" ) )
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" ) )

# **MAIN TEXT TABLES 6-10
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" & k.pub < 10 ) )
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" & k.pub < 10 ) )

# # **MAIN TEXT?
# # larger metas (interested in CI width here)
# make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" &
#                                                   k.pub >= 10 & k.pub <= 20))
# make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" &
#                                                   k.pub >= 10 & k.pub <= 20))


# ~ k=100  -------------------------------------------------

#**SUPP TABLES S1 - S4
# very large metas
make_both_winner_tables(.agg = agg %>% filter( Ytype == "cont-SMD" &
                                                 k.pub == 100))
make_both_winner_tables(.agg = agg %>% filter( Ytype == "bin-OR" &
                                                 k.pub == 100))

# sanity check: I think the reason Jeffreys has better MhatCoverNominal is only due to scens with expo population effects
#bm: not true??
make_both_winner_tables(.agg = agg %>% filter( Ytype == "cont-SMD" &
                                                 k.pub == 100 &
                                                 true.dist == "norm"))
make_both_winner_tables(.agg = agg %>% filter( Ytype == "bin-OR" &
                                                 k.pub == 100 &
                                                 true.dist == "norm"))

# # ~ Tau^2  -------------------------------------------------
# 
# # tau near zero
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", t2a == 0.0001))
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", t2a == 0.0001))
# 
# 
# # tau not near zero
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", t2a > 0.0001))
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", t2a > 0.0001))
# 
# #**SUPP TABLES S5 - S10
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", t2a <= 0.01))
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", t2a > 0.01))
# 
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", t2a <= 0.01))
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", t2a > 0.01))



# individual values
# make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.0025))
# make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.01))
# make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.04)) 



# ~ N.expr  -------------------------------------------------

make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", N.pretty %in% c("N ~ U(40, 400)", "N ~ U(2000, 4000)") ))
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", N.pretty %in% c("N = 40", "N = 400") ))

make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", N.pretty %in% c("N ~ U(40, 400)", "N ~ U(2000, 4000)") ))
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", N.pretty %in% c("N = 40", "N = 400") ))


# # individual values
# make_both_winner_tables(.agg = agg2 %>% filter(N.pretty == "N = 40"))
# make_both_winner_tables(.agg = agg2 %>% filter(N.pretty == "N = 400"))
# make_both_winner_tables(.agg = agg2 %>% filter( N.pretty == "N ~ U(40, 400)" ))
# make_both_winner_tables(.agg = agg2 %>% filter( N.pretty == "N ~ U(2000, 4000)" ))


# ~ Sensitivity analysis: When including scens with biased yi  -------------------------------------------------

# results are quite similar
setwd(data.dir)
list.files()
agge = fread( "agg_including_scens_biased_yi.csv") 
nuni(agge$scen.name)

agge2 = agge %>% filter(method %in% methods_for_table &
                          k.pub <= 20)



# overall winner tables
make_both_winner_tables(.agg = agge2 %>% filter( Ytype == "cont-SMD" ) )
make_both_winner_tables(.agg = agge2 %>% filter( Ytype == "bin-OR" ) )




# ONE-OFF PERFORMANCE STATS FOR PAPER -------------------------------------------------


### CI width comparisons

# SAVE?
# temp = agg2 %>% filter(method.pretty %in% methods_to_include) %>%
#   #filter(k.pub <= 5) %>%
#   filter(Ytype == "bin-OR") %>%
#   group_by(scen.name) %>%
#   mutate(jeff_wins1 = MhatCover[method.pretty == "Jeffreys2-central"] >= max(  MhatCover[method.pretty != "Jeffreys2-central"] ),
#          jeff_95 = MhatCover[method.pretty == "Jeffreys2-central"] >= .95,
#          jeff_wins2 = MhatWidth[method.pretty == "Jeffreys2-central"] <= min(  MhatWidth[method.pretty != "Jeffreys2-central"] ),
#          jeff_wins = jeff_wins1 * jeff_wins2)
# 
# mean(temp$jeff_95)
# mean(temp$jeff_wins1)
# mean(temp$jeff_wins2)
# mean(temp$jeff_wins)


# ~ Proportion of scens with nominal coverage  -------------------------------------------------



# Mhat
update_result_csv( name = paste( "Perc normal scens MhatCoverNominal Wald methods" ),
                   value = 100 * round( mean( agg2$MhatCoverNominal[ agg2$method.pretty %in% Wald_methods_pretty &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )


update_result_csv( name = paste( "Perc normal scens MhatCoverNominal MLE-profile" ),
                   value = 100 * round( mean( agg2$MhatCoverNominal[ agg2$method.pretty == "MLE-profile" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )


update_result_csv( name = paste( "Perc normal scens MhatCoverNominal Jeffreys1-central" ),
                   value = 100 * round( mean( agg2$MhatCoverNominal[ agg2$method.pretty == "Jeffreys1-central" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )


update_result_csv( name = paste( "Perc normal scens MhatCoverNominal Jeffreys2-central" ),
                   value = 100 * round( mean( agg2$MhatCoverNominal[ agg2$method.pretty == "Jeffreys2-central" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )



# Shat
( Q_methods = stringsWith(pattern =  "Qprofile", unique(agg2$method.pretty)) )
update_result_csv( name = paste( "Perc normal scens ShatCoverNominal Qprofile methods" ),
                   value = 100 * round( mean( agg2$ShatCoverNominal[ agg2$method.pretty %in% Q_methods &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )


update_result_csv( name = paste( "Perc normal scens ShatCoverNominal MLE-profile" ),
                   value = 100 * round( mean( agg2$ShatCoverNominal[ agg2$method.pretty == "MLE-profile" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )


update_result_csv( name = paste( "Perc normal scens ShatCoverNominal Jeffreys1-shortest" ),
                   value = 100 * round( mean( agg2$ShatCoverNominal[ agg2$method.pretty == "Jeffreys1-shortest" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )


update_result_csv( name = paste( "Perc normal scens ShatCoverNominal Jeffreys2-shortest" ),
                   value = 100 * round( mean( agg2$ShatCoverNominal[ agg2$method.pretty == "Jeffreys2-shortest" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )

update_result_csv( name = paste( "Perc normal scens ShatCoverNominal Jeffreys1-central" ),
                   value = 100 * round( mean( agg2$ShatCoverNominal[ agg2$method.pretty == "Jeffreys1-central" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )


update_result_csv( name = paste( "Perc normal scens ShatCoverNominal Jeffreys2-central" ),
                   value = 100 * round( mean( agg2$ShatCoverNominal[ agg2$method.pretty == "Jeffreys2-central" &
                                                                      agg2$true.dist == "norm" ] ), 2 ),
                   print = TRUE )



# ~ CI coverage and width comparisons: Bin-OR  -------------------------------------------------
# specifically compare to REML since all Wald-Qprofile methods were similar
x = CI_comparison(.agg = agg2 %>% 
                    filter(true.dist == "norm") %>%
                    filter(Ytype == "bin-OR"),
                  .target.method.pretty = "Jeffreys2-central",
                  .comparison.method.pretty = "REML-Wald-Qprofile")

update_result_csv( name = paste( "Bin-OR ", names(x) ),
                   value = as.numeric(x),
                   print = TRUE )

update_result_csv( name = "Sims - Mean perc narrower Jeffreys vs REML - Ybin",
                   value = perc_CI_narrower(.agg = agg2 %>% filter(Ytype == "bin-OR" & true.dist == "norm"),
                                            .target.method.pretty = "Jeffreys2-central",
                                            .comparison.method.pretty = "REML-Wald-Qprofile"),
                   print = TRUE )


update_result_csv( name = "Sims - Mean perc narrower Jeffreys vs REML - Ybin small metas",
                   value = perc_CI_narrower(.agg = agg2 %>% filter(Ytype == "bin-OR" & true.dist == "norm" & k.pub <= 5),
                                            .target.method.pretty = "Jeffreys2-central",
                                            .comparison.method.pretty = "REML-Wald-Qprofile"),
                   print = TRUE )


# ~ CI coverage and width comparisons: Cont-SMD  -------------------------------------------------
x = CI_comparison(.agg = agg2 %>% 
                    filter(k.pub > 5 & true.dist == "norm") %>%
                    filter(Ytype == "cont-SMD"),
                  .target.method.pretty = "Jeffreys2-central",
                  .comparison.method.pretty = "REML-Wald-Qprofile")

# **doesn't make sense to use Jeffreys2 for continuous outcomes since, in the k>5 scenarios where its coverage was always okay, it also doesn't improve efficiency
perc_CI_narrower(.agg = agg2 %>% filter(Ytype == "cont-SMD" & true.dist == "norm" & k.pub > 5),
                 .target.method.pretty = "Jeffreys2-central",
                 .comparison.method.pretty = "REML-Wald-Qprofile")






# SANITY CHECKS -------------------------

# ~ Compare to Kontopantelis results  -------------------------------------------------

make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" &
                                                  k.pub <= 5 &
                                                  t2a == 0.25 ) )


make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" &
                                                  k.pub <= 5 &
                                                  t2a == 0.0001 ) )

#***very interesting. Mhat coverage depends heavily on t2a for MLE-profile, but not for bayesmeta-joint-central! Both methods are fine for t2a = 0.0001, BUT only jeffreys is fine for t2a = 0.25.




# ~ Compare to Langan results  -------------------------------------------------

# for sanity checks, we include the scenarios in which even the individual studies were biased
# that only affects the inclusion of binary-Y scenarios, not continuous Y 

### try to reproduce Figure 1, upper left-hand panel
# this is one scen only
temp = agg %>% filter( Ytype == "cont-SMD" &
                         N.expr == "40",
                       k.pub == 2,
                       t2a == 0.04,
                       true.dist == "norm")
expect_equal( nuni(temp$scen.name), 1 )
make_both_winner_tables(.agg = temp,
                        display = "dataframe" )

# indeed, it is a negative bias, in contrast to Langan's observation of positive bias
temp$Shat
sqrt(temp$t2a)



### Is it true that PM has substantial positive bias per Langan?
temp = agg %>% filter(method == "PM")

summary(temp$ShatBias)
summary(temp$ShatBias / sqrt(temp$t2a))

# scens with >50% positive bias in Shat
temp2 = temp %>% filter( ShatBias > 0.5*sqrt(t2a) )
nrow(temp2)

# summarize scen params
CreateTableOne( dat = temp2,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2,
                strata = "Ytype" )

# across all scens
hist(temp$ShatBias / sqrt(temp$t2a))


# ~ Isolate an interesting scenario  -------------------------------------------------

x = agg2 %>% filter(t2a == 0 & ShatCover == 0)
unique(x$scen.name)

x2 = agg %>% filter(scen.name==1); View(x2)

# scen params only
x3 = x2[1, 1:12] 
library(constructive)
construct( as.data.frame(x3) )


# PLOTS -------------------------------------------------

# for the simple ggplotly versions
dp = agg2 %>% group_by(k.pub, t2a, method.pretty, Ytype.pretty, true.dist.pretty) %>%
  summarise_if(is.numeric, meanNA)

# includes k=100
dp100 = agg %>%
  #filter(k.pub == 100) %>%
  group_by(k.pub, t2a, method.pretty, Ytype.pretty, true.dist.pretty) %>%
  summarise_if(is.numeric, meanNA)

# ~ Plots by k and tau -------------------------------------------------


# ~ MhatMAE -------------------------------------------------


my_line_plot(.Yname = "MhatMAE",
             .agg = agg2,
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(mu) ~ " MAE") )',
             .jitter.width = 0.5)


# simple version for ggplotly
summary(dp$MhatMAE)
p = ggplot( data = dp, 
            aes( x = k.pub, 
                 y = MhatMAE,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  
  # slightly dodge line positions to avoid exact overlap:
  #geom_line( position=position_jitter(w=.5, h=0) ) + 
  geom_line() +
  scale_y_continuous(limits = c(0, .2)) +
  
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)

# ~ MhatRMSE -------------------------------------------------


my_line_plot(.Yname = "MhatRMSE",
             .agg = agg2,
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(mu) ~ " RMSE") )',
             .jitter.width = 0.5)


# simple version for ggplotly
summary(dp$MhatRMSE)
p = ggplot( data = dp, 
            aes( x = k.pub, 
                 y = MhatRMSE,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  
  # slightly dodge line positions to avoid exact overlap:
  #geom_line( position=position_jitter(w=.5, h=0) ) + 
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)



# ~ MhatCover -------------------------------------------------

my_line_plot(.Yname = "MhatCover",
             .agg = agg2,
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(mu) ~ " CI coverage") )',
             .jitter.width = 0)



# simple version for ggplotly
p = ggplot( data = dp, 
            aes( x = k.pub, 
                 y = MhatCover,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  geom_hline(yintercept = .95, lty = 1) +
  
  # slightly dodge line positions to avoid exact overlap:
  geom_line( position=position_jitter(w=.5, h=0) ) + 
  #geom_line() +
  
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)


# sanity check: understand coverage differences for k=100
# very interesting!
p = ggplot( data = dp100, 
            aes( x = k.pub, 
                 y = MhatCover,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  geom_hline(yintercept = .95, lty = 1) +
  
  # slightly dodge line positions to avoid exact overlap:
  geom_line( position=position_jitter(w=.5, h=0) ) + 
  #geom_line() +
  
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)


# ~ MhatWidth -------------------------------------------------

# zoom in to small k for legibility; all methods are basically the same for larger metas
my_line_plot(.Yname = "MhatWidth",
             #.agg = agg2 %>% filter(k.pub <= 10),
             .agg = agg2,
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(mu) ~ " CI width") )',
             .jitter.width = 0)

# simple version for ggplotly


p = ggplot( data = dp %>% filter(k.pub <= 10), 
            aes( x = k.pub, 
                 y = MhatWidth,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  
  # slightly dodge line positions to avoid exact overlap:
  #geom_line( position=position_jitter(w=.5, h=0) ) + 
  geom_line() +
  scale_y_log10() +
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)


# sanity check
agg2 %>% filter( Ytype == "cont-SMD" &
                   k.pub >= 10 & k.pub <= 20) %>%
  group_by(method.pretty.mu.inf) %>%
  summarise(meanNA(MhatWidth))

agg2 %>% filter( Ytype == "bin-OR" &
                   k.pub >= 10 & k.pub <= 20) %>%
  group_by(method.pretty.mu.inf) %>%
  summarise(meanNA(MhatWidth))

# SAVE - reproduce a single point (per method) on plot
agg2 %>% filter( Ytype == "bin-OR" &
                   k.pub == 10 &
                   t2a == 0.01) %>%
  group_by(method.pretty.mu.inf) %>%
  summarise(meanNA(MhatWidth))


# ~ ShatMAE -------------------------------------------------


my_line_plot(.Yname = "ShatMAE",
             .agg = agg2,
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(tau) ~ " MAE") )',
             .jitter.width = 0.5)


# simple version for ggplotly
p = ggplot( data = dp, 
            aes( x = k.pub, 
                 y = ShatMAE,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  
  # slightly dodge line positions to avoid exact overlap:
  #geom_line( position=position_jitter(w=.5, h=0) ) + 
  geom_line() +
  scale_y_continuous(limits = c(0, .4)) +
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)


# ~ ShatRMSE -------------------------------------------------


my_line_plot(.Yname = "ShatRMSE",
             .agg = agg2,
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(tau) ~ " RMSE") )',
             .jitter.width = 0.5)


# simple version for ggplotly
p = ggplot( data = dp, 
            aes( x = k.pub, 
                 y = ShatRMSE,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  
  # slightly dodge line positions to avoid exact overlap:
  #geom_line( position=position_jitter(w=.5, h=0) ) + 
  geom_line() +
  
  scale_y_continuous(limits = c(0, 0.4)) +
  
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)



# ~ ShatCover -------------------------------------------------

my_line_plot(.Yname = "ShatCover",
             .agg = agg2,
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(tau) ~ " CI coverage") )',
             .jitter.width = 0)



# simple version for ggplotly
p = ggplot( data = dp, 
            aes( x = k.pub, 
                 y = ShatCover,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  geom_hline(yintercept = .95, lty = 1) +
  
  # slightly dodge line positions to avoid exact overlap:
  #geom_line( position=position_jitter(w=.5, h=0) ) + 
  geom_line() +
  scale_y_continuous(limits = c(0.75, 1)) +
  
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)


# ~ ShatWidth -------------------------------------------------

#bm: rescale this one; it gets cut off at the bottom
# zoom in to small k for legibility; all methods are basically the same for larger metas
my_line_plot(.Yname = "ShatWidth",
             .agg = agg2 %>% filter(k.pub <= 10),
             .ggtitle = "",
             .ylab = 'bquote( bold( hat(tau) ~ " CI width") )',
             .jitter.width = 0)

# simple version for ggplotly
p = ggplot( data = dp %>% filter(k.pub <= 10), 
            aes( x = k.pub, 
                 y = ShatWidth,
                 color = method.pretty,
                 lty = method.pretty) ) + 
  
  # slightly dodge line positions to avoid exact overlap:
  #geom_line( position=position_jitter(w=.5, h=0) ) + 
  geom_line() +
  scale_y_log10() +
  facet_grid(t2a ~ true.dist.pretty + Ytype.pretty )

ggplotly(p)

# ~ Bias boxplots -------------------------------------------------

# ~~ MhatBias -----
my_boxplots(xName = "k.pub",
           yName = "MhatBias",
           hline = 0,
           xlab = "k",
           ylab = bquote( bold( hat(mu) ~ " bias") ),
           yTicks = round( seq(-0.04, 0.04, 0.01), 2 ),
           prefix = NA,
           
           write = TRUE,
           # by default, use all scenarios:
           .agg = agg2 )

# simple version for ggplotly with all methods:
quantile(agg2$MhatBias, c(0.05, 0.95), na.rm = TRUE)
p = ggplot( data = agg2 %>% filter(!is.na(method.pretty.est)),
            aes(x = k.pub,
                y = MhatBias,
                color = method.pretty.est,
                fill = method.pretty.est) ) +
  
  geom_boxplot(width=0.6,
               alpha = .5,
               outlier.shape = NA) +
  labs(color  = "Method", fill = "Method") +
  
  scale_y_continuous(limits = c(-0.04, 0.04)) +
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) +
  guides(color = guide_legend(nrow=2)) 
p




# ~~ ShatBias -------
my_boxplots(xName = "k.pub",
           yName = "ShatBias",
           hline = 0,
           xlab = "k",
           ylab = bquote( bold( hat(tau) ~ " bias") ),
           yTicks = round( seq(-0.4, 0.5, 0.02), 2 ),
           prefix = NA,
           
           write = TRUE,
           # by default, use all scenarios:
           .agg = agg2 )


# simple version for ggplotly with all methods:
quantile(agg2$ShatBias, c(0.05, 0.95), na.rm = TRUE)
p = ggplot( data = agg2 %>% filter(!is.na(method.pretty.est)),
            
            aes(x = as.factor(k.pub),
                y = ShatBias,
                color = method.pretty.est,
                fill = method.pretty.est) ) +
  
  geom_boxplot(width=0.6,
               alpha = .5,
               outlier.shape = NA) +
  labs(color  = "Method", fill = "Method") +
  
  scale_y_continuous(limits = c(-0.2, 0.6)) +
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) +
  guides(color = guide_legend(nrow=2)) +
  facet_grid( ~ Ytype.pretty )
p
ggplotly(p) %>% layout(boxmode = "group")


# ~~ Compare Jeffreys2 pmode, pmed, pmean  -------------------------------------------------

dp = agg2

dp$method.pretty = NA
dp$method.pretty[ dp$method == "jeffreys-pmode" ] = "Jeffreys2-mode"
dp$method.pretty[ dp$method == "jeffreys-pmean" ] = "Jeffreys2-mean"
dp$method.pretty[ dp$method == "jeffreys-pmed" ] = "Jeffreys2-median"

dp = dp %>% filter( !is.na(method.pretty) ) %>% droplevels()


### ShatBias
p = ggplot( data = dp,
            
            aes(x = as.factor(k.pub),
                y = ShatBias,
                color = method.pretty,
                fill = method.pretty) ) +
  
  geom_hline(yintercept = 0,
             lty = 2) +
  
  geom_boxplot(width=0.6,
               alpha = .5,
               outlier.shape = NA) +
  labs(color  = "Measure of central tendency", fill = "Measure of central tendency") +
  
  xlab("k") +
  ylab( bquote( bold( hat(tau) ~ " bias") ) ) +
  
  coord_cartesian(ylim = c(-0.5, 1.6)) +
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) +
  guides(color = guide_legend(nrow=2)) +
  facet_grid( ~ Ytype.pretty )
p

my_ggsave( name = "jeffreys2_ShatBias_point_estimates_comparison.pdf",
           .plot = p,
           .overleaf.dir = overleaf.dir.figs,
           .results.dir = results.dir,
           .width = 10,
           .height = 8)



### ShatRMSE
p = ggplot( data = dp,
            
            aes(x = as.factor(k.pub),
                y = ShatRMSE,
                color = method.pretty,
                fill = method.pretty) ) +
  
  geom_hline(yintercept = 0,
             lty = 2) +
  
  geom_boxplot(width=0.6,
               alpha = .5,
               outlier.shape = NA) +
  labs(color  = "Measure of central tendency", fill = "Measure of central tendency") +
  
  coord_cartesian(ylim = c(0, 1.6)) +
  
  xlab("k") +
  ylab( bquote( bold( hat(tau) ~ " RMSE") ) ) +
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) +
  guides(color = guide_legend(nrow=2)) +
  facet_grid( ~ Ytype.pretty )
p


my_ggsave( name = "jeffreys2_ShatRMSE_point_estimates_comparison.pdf",
           .plot = p,
           .overleaf.dir = overleaf.dir.figs,
           .results.dir = results.dir,
           .width = 10,
           .height = 8)


### ShatMAE
p = ggplot( data = dp,
            
            aes(x = as.factor(k.pub),
                y = ShatMAE,
                color = method.pretty,
                fill = method.pretty) ) +
  
  geom_hline(yintercept = 0,
             lty = 2) +
  
  # hide outliers:
  geom_boxplot(width=0.6,
               alpha = .5,
               outlier.shape = NA) +
  
  labs(color  = "Measure of central tendency", fill = "Measure of central tendency") +
  
  coord_cartesian(ylim = c(0, 1.5)) +
  
  xlab("k") +
  ylab( bquote( bold( hat(tau) ~ " MAE") ) ) +
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) +
  guides(color = guide_legend(nrow=2)) +
  facet_grid( ~ Ytype.pretty )
p


my_ggsave( name = "jeffreys2_ShatMAE_point_estimates_comparison.pdf",
           .plot = p,
           .overleaf.dir = overleaf.dir.figs,
           .results.dir = results.dir,
           .width = 10,
           .height = 8)


# SCEN 1384: Plots and stats (CI overcoverage despite better efficiency) -------------------------------------------------

# for binary Y, investigate the surprising finding that Jeffreys slightly over-covers, yet its CI is much narrower

# ~ Look for scens exhibiting this property  -------------------------------------------------
# wide form wrt methods:
wide_agg <- pivot_wider(agg2 %>% filter(method.pretty.mu.inf %in% methods_pretty_mu_inf),
                        names_from = method.pretty.mu.inf,
                        values_from = c(MhatCover, MhatWidth, MhatBias, MhatMAE,
                                        MLo, MHi),
                        names_sep = "_",
                        id_cols = all_of( c( "scen.name", param.vars.manip2 ) ) )

expect_equal( nrow(wide_agg), nuni(agg2$scen.name) )

# scens where Jeffreys over-covered but was narrower than REML:
scens = wide_agg$scen.name[ wide_agg$`MhatCover_Jeffreys2` > 0.95 & 
                              wide_agg$`MhatCover_REML-Wald` < 0.95 & 
                              wide_agg$`MhatWidth_Jeffreys2` < wide_agg$`MhatWidth_REML-Wald` ]

t = wide_agg %>% select( all_of( c( "scen.name", param.vars.manip2 ) ),
                         `MhatCover_Jeffreys2`, 
                         `MhatCover_REML-Wald`,
                         
                         `MhatWidth_Jeffreys2`, 
                         `MhatWidth_REML-Wald`,
                         
                         `MhatBias_Jeffreys2`, 
                         `MhatBias_REML-Wald`,
                         
                         `MhatMAE_Jeffreys2`, 
                         `MhatMAE_REML-Wald`,
                         
                         `MLo_Jeffreys2`, 
                         `MLo_REML-Wald`,
                         
                         `MHi_Jeffreys2`, 
                         `MHi_REML-Wald`) %>%
  filter(scen.name %in% scens) %>%
  arrange(`MhatCover_REML-Wald`) 
  #filter(scen.name == 1384)

if (use.View == TRUE) View(t)

# scen 1384 is striking; extract its scen params:
x = agg2 %>% filter(scen.name==1384)
x2 = x[1, 1:12] 
( scen_1384_params = constructive::construct( as.data.frame(x2) ) )

# I ran this scenario again locally for 100 sim reps to get iterate-level data

#$k=3$, binary $Y$,  $\mu = 0.5$, $\tau^2 = 0.04$,  normal population effects, $P(Y = 1 \mid X=0) = 0.05$,  and $N \sim U(2000, 4000)$
# temp: find a different scen
# this one was in the manuscript previously
temp = agg %>% filter(k.pub == 3 &
                 Mu == 0.5 &
                 t2a == 0.04 &
                 true.dist == "norm" &
                 p0 == 0.05 &
                 Ytype == "bin-OR" &
                 N.pretty == "N ~ U(2000, 4000)")

nrow(temp)
temp$scen.name


# ~ Look at individual iterates for scen 1384  -------------------------------------------------

# scen 1384 data
setwd(data.dir)
s2 = fread("pretty_long_results_job_1384.csv")


### MhatBias vs. MhatWidth
# **super interesting!!

mean(is.na(s2$MhatCover))

# plotting df
methods_for_plot = rev( c("DL-Wald", "REML-Wald", "Exact", "Jeffreys2", "Jeffreys1") )

s2p = s2 %>% filter( !is.na(MhatCover) &
                       !is.na(method.pretty.mu.inf) &
                       method.pretty.mu.inf %in% methods_for_plot ) %>%
  droplevels()
  
# reorder methods
correct.order = rev( c("DL-Wald", "REML-Wald", "Exact", "Jeffreys1", "Jeffreys2") )
s2p$method.pretty.mu.inf = factor(s2p$method.pretty.mu.inf, levels = correct.order)
any(is.na(s2p$method.pretty.mu.inf))

# find good x-axis limits
summary( s2p$MhatBias )
quantile(s2p$MhatBias, 0.95, na.rm = TRUE)
xmin = -0.10
xmax = 0.20

# find good y-axis limits
summary( s2p$MhatWidth )
ymin = 0
ymax = 5

# same colors as in analyze_sims_helper.R for prettiness
.colors = rev( c("#246105",
                 "black",
                 "#CC9808",
                 "#F2340E",
                 "#E075DB" ) )


p = ggplot( data = s2p, 
            aes( x = MhatBias,
                 y = MhatWidth,
                 color = as.factor(MhatCover) ) ) +
  geom_point(alpha = 0.5) +
  facet_wrap( ~ method.pretty.mu.inf,
              nrow = 2) +
  
  scale_color_manual( values = c("red", "black") ) +
  
  scale_x_continuous(limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, 0.05) ) +
  
  scale_y_continuous( limits = c(ymin, ymax),
                      breaks = seq(ymin, ymax, .5) ) +
  
  labs(color  = bquote(CI ~ includes ~ mu) ) +
  xlab( bquote(hat(mu) ~ bias) ) +
  ylab( "CI width" ) +
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) 
p




my_ggsave( name = "scen_1384_bias_vs_CI_width.pdf",
           .plot = p,
           .overleaf.dir = overleaf.dir.figs,
           .results.dir = results.dir,
           .width = 10,
           .height = 8)



### CI asymmetry across scenarios
p = ggplot( data = s2 %>% filter(method.pretty == "Jeffreys"), 
            aes( x = CI_asy) ) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_density(size = 1.2) +
  
  theme_bw(base_size = 16) +
  
  xlab( "CI asymmetry" ) +
  ylab("Density") +
  scale_y_continuous(breaks = NULL) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         axis.ticks.y = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) 

my_ggsave( name = "scen_1384_CI_asymmetry.pdf",
           .plot = p,
           .overleaf.dir = overleaf.dir.figs,
           .results.dir = results.dir,
           .width = 8,
           .height = 6)


### One-off stats about this scenario

# parameters, hard-coded in Supplement
scen_1384_params


temp = agg2 %>% filter(scen.name == 1384 & method.pretty.mu.inf %in% methods_pretty_mu_inf)


update_result_csv( name = "Scen 1384 Jeffreys2 MhatCover",
                   value = round( 100 * as.numeric( temp %>% filter(method.pretty.mu.inf == "Jeffreys2") %>%
                                                      select(MhatCover) ), 6 ),
                   print = TRUE )

# all Wald
coverages = round( 100 * temp$MhatCover[ temp$method.pretty.mu.inf %in% Wald_methods_pretty ] )
expect_equal( nuni(coverages), 1 )  # otherwise doesn't make to summarize by a single number as below
update_result_csv( name = "Scen 1384 Wald methods MhatCover",
                   value = unique(coverages),
print = TRUE )


update_result_csv( name = "Scen 1384 Jeffreys2 MhatWidth",
                   value = round( as.numeric( temp %>% filter(method.pretty.mu.inf == "Jeffreys2") %>%
                                                select(MhatWidth) ), 2 ),
                   print = TRUE )


# all other methods
round( temp$MhatWidth[ temp$method.pretty.mu.inf %in% Wald_methods_pretty ], 2 )


# ~ Line plots of multiple outcomes (not in use)  -------------------------------------------------

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

if (use.View = TRUE) View(t)



# LOOK FOR SCENS WITH MORE VARIABILITY ON PERFORMANCE METRICS  -------------------------------------------------

# scens with meaningful differences across methods in the outcomes

methods_for_table = c(
  "ML", "REML", "DL", "PM", "DL2", "exact","MLE-profile",
  "bayesmeta-tau-central", 
  "bayesmeta-tau-shortest",
  "bayesmeta-joint-shortest",
  "bayesmeta-joint-central")

# alternate: range instead of SD
t = agg2 %>% group_by(scen.name) %>%
  # important to avoid bas methods like pmean, pmed, etc.
  filter(method %in% methods_for_table) %>%
  summarise( MhatBiasRange = diff( range(MhatBias, na.rm = TRUE) ),
             MhatMAERange = diff( range(MhatMAE, na.rm = TRUE) ),
             MhatCoverRange = diff( range(MhatCover, na.rm = TRUE) ),
             MhatRMSERange = diff( range(MhatRMSE, na.rm = TRUE) ),
             
             ShatBiasRange = diff( range(ShatBias, na.rm = TRUE) ),
             ShatMAERange = diff( range(ShatMAE, na.rm = TRUE) ),
             ShatCoverRange = diff( range(ShatCover, na.rm = TRUE) ),
             ShatRMSERange = diff( range(ShatRMSE, na.rm = TRUE) ) )


#**performance metrics for paper
update_result_csv( name = "Max MhatRMSERange across scens",
                   value = max(t$MhatRMSERange),
                   print = TRUE )
update_result_csv( name = "Max MhatMAERange across scens",
                   value = max(t$MhatMAERange),
                   print = TRUE )
update_result_csv( name = "Max MhatBiasRange across scens",
                   value = max(t$MhatBiasRange),
                   print = TRUE )



# *much more variability on Shat performance metrics than on Mu
summary(t$ShatMAERange)
summary(t$ShatCoverRange)
summary(t$ShatRMSERange)
summary(t$ShatBiasRange)


# mark scens as important if methods vary on any of these characteristics
t$scen_important_MhatCoverRange = t$MhatCoverRange > 0.1
t$scen_important_ShatCoverRange = t$ShatCoverRange > 0.1

t$scen_important_MhatRMSERange = t$MhatRMSERange > 0.25
t$scen_important_ShatRMSERange = t$ShatRMSERange > 0.25

t$scen_important_MhatMAERange = t$MhatMAERange > 0.1
t$scen_important_ShatMAERange = t$ShatMAERange > 0.1

# percent of scens that were important by each metric
as.data.frame( t %>% select( namesWith("scen_important", t) ) %>%
  summarise_all(meanNA) )


# # **add the importance vars to agg
# agg = left_join(x = agg,
#                 y = t %>% select(scen.name, scen_important),
#                 by = "scen.name")

# add these vars to agg
agg = left_join(x = agg,
                y = t %>% select(scen.name, namesWith("Range", t) ),
                by = "scen.name")


# EXPLORE PREDICTORS OF HIGHER RANGE ON EACH OUTCOME VAR  -------------------------------------------------

# must run code in previous section to get the range variables

regressions.from.scratch = TRUE

# remove p0 because causes error "contrasts can be applied only to factors with 2 or more levels" in lm
#  because p0 doesn't vary conditional on Ytype = cont
param.vars.manip3 = drop_vec_elements(param.vars.manip2, "p0")
param.vars.manip3[ param.vars.manip3 == "N.expr" ] = "N.pretty"

( range_names = namesWith("scen_important", agg) )

( t1 = performance_regressions(.agg = agg,
                               Ynames = range_names,
                               covariates = param.vars.manip3 ) )
View(t1)

