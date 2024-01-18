# NOTES ----------------------------------------------------


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
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/JTE (Jeffreys tau estimation) Overleaf/R_objects/figures"
# test it
setwd(overleaf.dir.figs)


setwd(code.dir)
source("analyze_sims_helper_JTE.R")
source("helper_JTE.R")  # for lprior(), etc.


# ** LINE PLOT: TAU FOR DIFFERENT N.EXPR  -------------------------------------------------

# Idea: panels are different values of k; lines are different sei distributions taken from the sims

# as in the simulation study

k.pub = c(2, 5, 10, 100)


# list of all results (one for each k.pub)
resl = list()
# list of plots (one for each k.pub)
pl = list()

for ( i in 1:length(k.pub) ) {
  res = prior_plot_one_k(.k = k.pub[i])
  resl[[i]] = res
  pl[[i]] = res$plot
}

# prior shape is really almost the same regardless of k
pl[[1]]
pl[[2]]
pl[[3]]
pl[[4]]

# save just one of them
plot = pl[[3]] + ggtitle("")

my_ggsave(name = "prior_plot.pdf",
          .plot = plot,
          .width = 10,
          .height = 6,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)


# # CONTOUR PLOT - not super useful -------------------------------------------------
# 
# # calculate the prior for different values
# dp = expand_grid( .mu = seq(-1, 1, 0.01),
#                   .tau = seq(0, 1, 0.05) )
# nrow(dp)
# 
# 
# 
# dp = dp %>% rowwise() %>%
#   mutate( prior.val = get_lprior(mu = .mu, tau = .tau, sei = 1)  )
# #bm
# 
# # set up colors for contours
# get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
# myColors = get_colors(n=15)  # chose 11 based on errors from ggplot if it was fewer
# 
# ### Contour plot ###
# p1 = ggplot( data = dp, 
#              aes(x = .mu,
#                  y = .tau,
#                  z = prior.val) ) +
#   
#   geom_contour_filled() +
#   
#   # close, but not enough colors
#   scale_fill_manual(values = myColors) +
#   
#   geom_contour(color = "white") +
#   
#   xlab( bquote(mu) ) +
#   ylab( bquote(tau) ) +
#   
#   geom_vline( xintercept = 0, lty = 2 ) +
#   
#   scale_y_continuous(breaks = seq( min(dp$.tau), max(dp$.tau), 0.25),
#                      limits = c( min(dp$.tau), max(dp$.tau) ) ) +
#   
#   theme_bw(base_size = 16) +
#   theme(text = element_text(face = "bold"),
#         axis.title = element_text(size=20),
#         legend.position = "none")
# 
# p1
# 
# 
# ### Line plot: continuous mu, discrete tau
# 
# dp3 = expand_grid( .mu = seq(-1, 1, 0.05),
#                    .tau = c(0, 0.25, .5, 0.75) )
# nrow(dp3)
# dp3 = dp3 %>% rowwise() %>%
#   mutate( prior.val = get_lprior(mu = .mu, tau = .tau, sei = 1)  )
# 
# 
# p3 = ggplot( data = dp3, 
#              aes(x = .mu,
#                  y = prior.val,
#                  color = as.factor(.tau) ) ) +
#   
#   geom_line() +
#   
#   xlab( bquote(.mu) ) +
#   ylab( "Log prior" ) +
#   
#   geom_vline( xintercept = 0, lty = 2 ) +
#   
#   # scale_y_continuous(breaks = seq( min(dp$.tau), max(dp$.tau), 0.25),
#   #                    limits = c( min(dp$.tau), max(dp$.tau) ) ) +
#   
#   theme_bw(base_size = 16) +
#   theme(text = element_text(face = "bold"),
#         axis.title = element_text(size=20) )
# 
# p3
# 
