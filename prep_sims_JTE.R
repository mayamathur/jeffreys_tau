
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
# renv::snapshot()


# ~~ User-specified global vars -------------------------
# no sci notation
options(scipen=999)

stitch.from.scratch = FALSE


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


setwd(code.dir)
source("helper_JTE.R")
source("analyze_sims_helper_JTE.R")




# READ IN AGGREGATED DATA FROM CLUSTER -------------------------------------------------


# ~ Basic prep -------------------------------------------------
setwd(data.dir)
aggo = fread("aggo.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("aggo.csv")$mtime
nrow(aggo) / nuni(aggo$method)  # number of scens that are done; 2496 if sims are done

# add fancy variables for plotting, etc.
agg = wrangle_agg_local(aggo)

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names(.agg = agg)


# ~ Individual studies should be unbiased -------------------------------------------------

# look for scens where even the individual studies are biased for Mu
#  e.g., because of very rare binary Y with small N
summary( abs(agg$sancheck_mean_yi - agg$Mu) )

# flag scens where yi has bias > 0.05
agg$exclude_scen_biased_yi = abs(agg$sancheck_mean_yi - agg$Mu) > 0.05
mean(agg$exclude_scen_biased_yi)  # percent of scens

# summarize scen params for these ones
# not surprisingly, the bad scens are exclusively binary Y, and mostly ones with N=40 or N ~ Unif(40,400)
#  though spread across a variety of p0 values
agg_bad = agg %>% filter(exclude_scen_biased_yi == TRUE)
CreateTableOne( dat = agg_bad,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2 )

table(agg_bad$N.expr )



#***exclude these scens going forward
agg = agg %>% filter(exclude_scen_biased_yi == FALSE)
CreateTableOne( dat = agg,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2,
                strata = "Ytype")


# Write prepped datasets -------------------------------------------------

setwd(data.dir)
fwrite(agg, "agg.csv")
fwrite(agg_bad, "agg_excluded_scens_biased_yi.csv")


# PREP ITERATE-LEVEL DATA FOR SCEN 1072  -------------------------------------------------

# this will be a little slow (1-2 min)
setwd(data.dir)
s = fread("stitched.csv")

# similar prep work to wrangle_agg_local
s2 = s %>% filter(scen.name == 1072 &
                    method %in% c("REML", "DL", "PM", "DL2", "robu", "jeffreys-pmode") )


s2 = s2 %>% mutate( CI_asy = (MHi - Mhat) / (Mhat - MLo),
                    MhatBias = Mhat - Mu,
                    MhatWidth = MHi - MLo,
                    MhatCover = (MHi >= Mu & MLo <= Mu) )

s2$method.pretty = s2$method 
s2$method.pretty[ s2$method == "robu" ] = "RVE"
s2$method.pretty[ s2$method == "jeffreys-pmode" ] = "Jeffreys"
table(s2$method, s2$method.pretty)

correct.order = c( "Jeffreys", "DL", "DL2", "REML", "PM", "RVE")
s2$method.pretty = factor(s2$method.pretty, levels = correct.order)
levels(s2$method.pretty)


fwrite( s2, "stitched_scen_1072.csv" )




# NOT IN USE (FROM OTHER PROJECTS):
# MERGE 4 SIMULATION DATASETS (ITERATE LEVEL) -------------------------------------------------

# # bind the stitched files
# for (i in 1:length(data.dir.suffixes) ) {
#   
#   .dir = data.dir.suffixes[i]
#   
#   setwd(data.dir)
#   setwd(.dir)
#   
#   s.chunk = fread("stitched.csv")
#   
#   summary(s.chunk$scen.name)
#   
#   if (i == 1) {
#     s = s.chunk
#     
#     # just for sanity checks
#     sanity = data.frame(sim.env = s.chunk$sim.env[1],
#                         methods = s.chunk$rep.methods[1],
#                         n.methods = nuni(s.chunk$method),
#                         n.scens = nuni(s.chunk$scen.name))
#     
#   } else {
#     # *need to bind_rows here to fill in NA columns (e.g., vars that don't apply for stefan sim env)
#     # hence approach of directly binding the iterate-level data before aggregating
#     s = bind_rows(s, s.chunk)
#     
#     sanity = bind_rows(sanity, 
#                        data.frame(sim.env = s.chunk$sim.env[1],
#                                   methods = s.chunk$rep.methods[1],
#                                   n.methods = nuni(s.chunk$method),
#                                   n.scens = nuni(s.chunk$scen.name)) )
#   }
#   
#   setwd(results.dir)
#   fwrite(s, "stitched_merged.csv")
#   fwrite(sanity, "sanity.csv")
#   
# } # end loop over data.dir.suffixes




# # STITCH FROM SCRATCH -------------------------------------------------
# 
# if ( stitch.from.scratch == TRUE ) {
#   setwd(data.dir)
#   s = fread("stitched.csv")
#   
#   aggo = make_agg_data(s,
#                        expected.sim.reps = 1000) 
#   
#   aggo = make_agg_data(s[1:20000,],
#                        expected.sim.reps = 1000) 
#   
#   
#   # sanity check:
#   # should have 1 row per scen-method combo
#   n.scens = 300
#   n.methods = nuni(aggo$method)
#   expect_equal( sum( n.scens * n.methods ),
#                 nrow(aggo) )
#   setwd(data.dir); fwrite(aggo, "agg.csv")
#   
#   
#   # add fancy variables for plotting, etc.
#   agg = wrangle_agg_local(aggo)
#   setwd(data.dir); fwrite(agg, "agg.csv")
#   
# }



