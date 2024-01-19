
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




# STITCH FROM SCRATCH -------------------------------------------------

if ( stitch.from.scratch == TRUE ) {
  setwd(data.dir)
  s = fread("stitched.csv")
  
  aggo = make_agg_data(s,
                       expected.sim.reps = 1000) 
  
  aggo = make_agg_data(s[1:20000,],
                       expected.sim.reps = 1000) 
  
  
  # sanity check:
  # should have 1 row per scen-method combo
  n.scens = 300
  n.methods = nuni(aggo$method)
  expect_equal( sum( n.scens * n.methods ),
                nrow(aggo) )
  setwd(data.dir); fwrite(aggo, "agg.csv")
  
  
  # add fancy variables for plotting, etc.
  agg = wrangle_agg_local(aggo)
  setwd(data.dir); fwrite(agg, "agg.csv")
  
}


# ALTERNATIVE: READ IN AGG DATA FROM CLUSTER -------------------------------------------------

setwd(data.dir)
aggo = fread("aggo.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("aggo.csv")$mtime
nrow(aggo) / nuni(aggo$method)  # number of scens that are done; 2496 if sims are done

# add fancy variables for plotting, etc.
agg = wrangle_agg_local(aggo)
setwd(data.dir); fwrite(agg, "agg.csv")




