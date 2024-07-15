
# PRELIMINARIES ----------------------------------------------------

# rm(list=ls())

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



setwd(code.dir)
source("analyze_sims_helper_JTE.R")
source("helper_JTE.R")  # for lprior(), etc.




# CHECK SMD VARIANCE SIMPLIFICATION  -------------------------------------------------

# sanity check:


# from debug(escalc):
# For measure="SMD", one can choose between vtype="LS" (the default) for the usual large-sample approximation to compute the sampling variances (equation 8 in Hedges, 1982),
# if (vtype[i] == "LS") 
#   vi[i] <- 1/n1i[i] + 1/n2i[i] + yi[i]^2/(2 * 
#                                             npi[i])
# where npi = total N

N = 10
Mu = -4

d = sim_meta( 
  Mu = Mu,
  t2a = 0,
  true.dist = "norm",
  
  N.expr = N,
  Ytype = "cont-SMD",
  p0 = NA,
  
  k.pub = 500)


# manual approximation of escalc formula above for case of E[N1] = E[N2]
my_sei = sqrt( (8+Mu^2)/(2*N) )

summary(d$sei - my_sei)


(8+Mu^2)/(2*N)  # manual simplification of escalc formula above for case of E[N1] = E[N2]
# yes, matches :)



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
