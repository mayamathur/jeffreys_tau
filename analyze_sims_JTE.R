

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

# should all the "View()" things be called?
# they open new RStudio tabs, so can be useful to skip these sanity checks
use.View = FALSE

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




# ~~ Get agg data -------------------------

# if only analyzing a single sim environment (no merging):
setwd(data.dir)
# "agg all" because we will later exclude scenarios whose parameters led to within-study bias
agga = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agga) / nuni(agga$method)


# drop any "NA" methods (i.e., ones that didn't get labeled in wrangle_agg_local)
agga = agga %>% filter( method.pretty != "" )
table(agga$method.pretty)

# look at number of actual sim reps
table(agga$sim.reps.actual)



# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names(.agg = agga)


# summarize scen params
CreateTableOne( dat = agga,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2,
                strata = "Ytype" )


# ~~ Check runtimes of sbatch files -------------------------

# mean runtimes within scenarios - HOURS
summary(agga$doParallelSeconds/60^2) 

# 95th quantile of runtime within scens - HOURS
summary(agga$doParallelSecondsQ95/60^2) 



                        

# SANITY CHECKS ON DATA GENERATION -------------------------


# ~ Individual studies should be unbiased  -------------------------------------------------

# look for scens where even the individual studies are biased for Mu
#  e.g., because of very rare binary Y with small N
summary( abs(agga$sancheck_mean_yi - agga$Mu) )

ind = which( abs(agga$sancheck_mean_yi - agga$Mu) > 0.1 )
length(ind)/nrow(agga)  # percent of scens

# summarize scen params for these ones
CreateTableOne( dat = agga[ind,],
                vars = param.vars.manip2,
                factorVars = param.vars.manip2)

table(agga[ind,"N.expr"] )
# not surprisingly, the bad scens are exclusively binary Y, and almost exclusively ones with N=40
#  though spread across a variety of p0 values


# now with a more stringent threshold for bias in yi
ind = which( abs(agga$sancheck_mean_yi - agga$Mu) > 0.05 )
ind
length(ind)/nrow(agga)  # percent of scens

# summarize scen params for these ones
CreateTableOne( dat = agga[ind,],
                vars = param.vars.manip2,
                factorVars = param.vars.manip2 )

table(agga[ind,"N.expr"] )
# not surprisingly, the bad scens are exclusively binary Y, and almost exclusively ones with N=40
#  though spread across a variety of p0 values


#***exclude these scens going forward
agg = agga[-ind,]
CreateTableOne( dat = agg,
                vars = param.vars.manip2,
                factorVars = param.vars.manip2,
                strata = "Ytype")

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

if (use.View = TRUE) View(t)




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


# ******** WINNER TABLES -------------------------

# create the base dataset from which to filter all winner tables
#agg2 = agg %>% filter( scen_important == TRUE & Ytype == .Ytype )
agg2 = agg %>% filter( k.pub <= 20)


dim(agg2)
# summarize scen params
CreateTableOne( dat = agg2[ !duplicated(agg2$scen.name) ],
                vars = param.vars.manip2,
                factorVars = param.vars.manip2 )


# ~ Overall  -------------------------------------------------

# **MAIN TEXT TABLES 2-5
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" ) )
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" ) )

# **MAIN TEXT TABLES 6-10
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "cont-SMD" & k.pub <= 5 ) )
make_both_winner_tables(.agg = agg2 %>% filter( Ytype == "bin-OR" & k.pub <= 5 ) )



# ~ k  -------------------------------------------------

#**SUPP TABLES S1 - S4
# very small metas
make_both_winner_tables(.agg = agg %>% filter( Ytype == "cont-SMD", k.pub == 100))
make_both_winner_tables(.agg = agg %>% filter( Ytype == "bin-OR", k.pub == 100))


# ~ Tau^2  -------------------------------------------------

#**SUPP TABLES S5 - S10
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", t2a <= 0.01))
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", t2a > 0.01))

make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", t2a <= 0.01))
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", t2a > 0.01))



# individual values
# make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.0025))
# make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.01))
# make_both_winner_tables(.agg = agg2 %>% filter(t2a == 0.04)) 


# ~ True effect distribution  -------------------------------------------------
# makes little difference for all methods
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", true.dist == "norm"))
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", true.dist == "expo"))

make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", true.dist == "norm"))
make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", true.dist == "expo"))

# # what about high-t2a scens? that's where I'd expect dist to matter more
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", true.dist == "norm", t2a == 0.25))
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "cont-SMD", true.dist == "expo", t2a == 0.25))
# 
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", true.dist == "norm", t2a == 0.25))
# make_both_winner_tables(.agg = agg2 %>% filter(Ytype == "bin-OR", true.dist == "expo", t2a == 0.25))

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



# ONE-OFF PERFORMANCE STATS FOR PAPER -------------------------------------------------


### CI width comparisons

# these use only the k <= 20 scens


#@UNDERLYING FN NEEDS SANITY CHECKS
update_result_csv( name = "Sims - Mean perc narrower Jeffreys vs winning other method - Ycont",
                   value = perc_CI_narrower(.agg = agg2 %>% filter(Ytype == "cont-SMD")),
                   print = TRUE )
update_result_csv( name = "Sims - Mean perc narrower Jeffreys vs winning other method - Ycont small metas",
                   value = perc_CI_narrower(.agg = agg2 %>% filter(k.pub <= 5 & Ytype == "cont-SMD")),
                   print = TRUE )



update_result_csv( name = "Sims - Mean perc narrower Jeffreys vs winning other method - Ybin",
                   value = perc_CI_narrower(.agg = agg2 %>% filter(Ytype == "bin-OR")),
                   print = TRUE )

update_result_csv( name = "Sims - Mean perc narrower Jeffreys vs winning other method - Ybin small metas",
                   value = perc_CI_narrower(.agg = agg2 %>% filter(k.pub <= 5 & Ytype == "bin-OR")),
                   print = TRUE )


### How asymmetric is Jeffreys CI for binary outcomes?

# for binary Y, investigate the surprising finding that Jeffreys slightly over-covers, yet its CI is much narrower
# wide form wrt methods:
w = pivot_wider(agg, names_from = method, values_from = MhatCover)


wide_agg <- pivot_wider(agg,
                        names_from = method,
                        values_from = c(MhatCover, MhatWidth, MhatBias, MhatAbsBias,
                                        MLo, MHi),
                        names_sep = "_",
                        id_cols = all_of( c( "scen.name", param.vars.manip2 ) ) )

expect_equal( nrow(wide_agg), nuni(agg$scen.name) )
# View( wide_agg %>% select(`MhatCover_jeffreys-pmode`, `MhatCover_REML`) )

# scens where Jeffreys over-covered but was narrower than REML:
scens = wide_agg$scen.name[ wide_agg$`MhatCover_jeffreys-pmode` > 0.95 & 
                              wide_agg$`MhatCover_REML` < 0.95 & 
                              wide_agg$`MhatWidth_jeffreys-pmode` < wide_agg$`MhatWidth_REML` ]

t = wide_agg %>% select( all_of( c( "scen.name", param.vars.manip2 ) ),
                         `MhatCover_jeffreys-pmode`, 
                         `MhatCover_REML`,
                         
                         `MhatWidth_jeffreys-pmode`, 
                         `MhatWidth_REML`,
                         
                         `MhatBias_jeffreys-pmode`, 
                         `MhatBias_REML`,
                         
                         `MhatAbsBias_jeffreys-pmode`, 
                         `MhatAbsBias_REML`,
                         
                         `MLo_jeffreys-pmode`, 
                         `MLo_REML`,
                         
                         `MHi_jeffreys-pmode`, 
                         `MHi_REML`) %>%
  filter(scen.name %in% scens) 
#filter(scen.name == 1072)

if (use.View = TRUE) View(t)
# scen 1072 is striking

# asymmetry
t = agg %>% mutate(CI_asy = (MHi - Mhat) / (Mhat - MLo) )

temp = t %>% group_by(method.pretty) %>%
  filter(method.pretty %in% methods.to.show & scen.name == 1072) %>%
  summarise( mean(CI_asy))

temp

### look at individual iterates for scen 1072

# commented out because slow
#setwd(data.dir)
#s = fread("stitched.csv")
#s2 = s %>% filter(scen.name == 1072)
#fwrite( s2, "stitched_scen_1072.csv" )

# prep the data
setwd(data.dir)
s2 = fread("stitched_scen_1072.csv") %>% filter(method %in% c("REML", "DL", "PM", "DL2", "robu", "jeffreys-pmode"))

s2 = s2 %>% mutate( CI_asy = (MHi - Mhat) / (Mhat - MLo),
                    MhatBias = Mhat - Mu,
                    MhatWidth = MHi - MLo,
                    MhatCover = (MHi >= Mu & MLo <= Mu) )

s2$method.pretty = s2$method # temporarily not relabeling them
s2$method.pretty[ s2$method == "robu" ] = "RVE"
s2$method.pretty[ s2$method == "jeffreys-pmode" ] = "Jeffreys"
table(s2$method, s2$method.pretty)

correct.order = c( "Jeffreys", "DL", "DL2", "REML", "PM", "RVE")
s2$method.pretty = factor(s2$method.pretty, levels = correct.order)
levels(s2$method.pretty)


#bm
# **super interesting!!
p = ggplot( data = s2, 
        aes( x = MhatBias,
             y = MhatWidth,
             color = MhatCover) ) +
  geom_point(alpha = 0.5) +
  facet_wrap( ~ method.pretty,
              nrow = 2) +
  
  scale_color_manual( values = c("red", "black") ) +
  scale_y_continuous( limits = c(0, 2.6),
                      breaks = seq(0, 3, .5) ) +
  
  labs(color  = bquote(CI ~ includes ~ mu) ) +
  xlab( bquote(hat(mu) ~ bias) ) +
  ylab( "CI width" ) +
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" ) 

my_ggsave( name = "scen_1072_bias_vs_CI_width.pdf",
           .plot = p,
           .overleaf.dir = overleaf.dir.figs,
           .results.dir = results.dir,
           .width = 10,
           .height = 8)



### 
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

my_ggsave( name = "scen_1072_CI_asymmetry.pdf",
           .plot = p,
           .overleaf.dir = overleaf.dir.figs,
           .results.dir = results.dir,
           .width = 8,
           .height = 6)

#bm: save these gorgeous plots! :D



# SANITY CHECKS -------------------------

# ~ Investigate binary-Y CI -------------------------------------------------

# why is it that Jeffreys over-covers for binary Y, yet is much narrower?

# ~ Compare to Langan results  -------------------------------------------------

# for sanity checks, we include the scenarios in which even the individual studies were biased
# that only affects the inclusion of binary-Y scenarios, not continuous Y 

### overall
# indeed, all methods have substantial positive bias
make_both_winner_tables(.agg = agga %>% filter( Ytype == "cont-SMD" ),
                        display = "dataframe")
# c.f. excluding the problematic scen combos:
make_both_winner_tables(.agg = agg %>% filter( Ytype == "cont-SMD" ),
                        display = "dataframe" )


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



### is it true that PM has substantial positive bias per Langan?
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

x = agg2 %>% filter(t2a == 0.0001 & ShatCover == 0)
unique(x$scen.name)

x2 = agg %>% filter(scen.name==105); View(x2)

# scen params only
x3 = x2[1, 1:11] 
library(constructive)
construct( as.data.frame(x3) )


# PLOTS -------------------------------------------------


# ~ Boxplots for (signed) bias -------------------------------------------------

aggp = agg2 %>% filter(method.pretty %in% methods.to.show )
unique(aggp$method.pretty)

# reorder methods
correct.order = c( "Jeffreys", "DL", "DL2", "REML", "PM", "RVE")
aggp$method.pretty = factor(aggp$method.pretty, levels = correct.order)
levels(aggp$method.pretty)

aggp$Ytype.pretty = NA
aggp$Ytype.pretty[ aggp$Ytype == "cont-SMD" ] = "Continuous Y"
aggp$Ytype.pretty[ aggp$Ytype == "bin-OR" ] = "Binary Y"
aggp$Ytype.pretty = factor(aggp$Ytype.pretty, levels = c("Continuous Y", "Binary Y") )

# same colors as in applied example for prettiness
.colors = c("#F2340E",
            "#0E96F0",
            "#845699",
            "#0F5A8C",
            "#6D9956",
            "#D18350")


my_violins(xName = "k.pub",
           yName = "MhatBias",
           hline = 0,
           xlab = "k",
           ylab = bquote( bold( hat(mu) ~ " bias") ),
           yTicks = round( seq(-0.04, 0.04, 0.01), 2 ),
           prefix = NA,
           colors = .colors,
           
           write = TRUE,
           # by default, use all scenarios:
           .agg = aggp )

my_violins(xName = "k.pub",
           yName = "MhatWidth",
           hline = 0,
           xlab = "k",
           ylab = bquote( bold( hat(mu) ~ " CI width") ),
           #yTicks = round( seq(-0.04, 0.04, 0.01), 2 ),
           prefix = NA,       
           colors = .colors,
           write = TRUE,
           # by default, use all scenarios:
           .agg = aggp )





my_violins(xName = "k.pub",
           yName = "ShatBias",
           hline = 0,
           xlab = "k",
           ylab = bquote( bold( hat(tau) ~ " bias") ),
           yTicks = round( seq(-0.40, 0.60, 0.1), 2 ),
           prefix = NA,
           colors = .colors,
           write = TRUE,
           # by default, use all scenarios:
           .agg = aggp )








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

