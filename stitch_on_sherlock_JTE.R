


# run this interactively in ml load R or via:
#   sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/job_stitch.sbatch
# scontrol show job 33834701
# look at its out file:
# cd /home/groups/manishad/SAPH
# cd /home/users/mmathur
# less rm_stitch.out

# for non-huge simulations, can often run this script interactively in a higher-memory
#  Sherlock session:
# ml load R/4.1.2
# srun --mem=32G --time=3:00:00 --pty bash
# R


# to be run by stitch.sbatch or manually
# To quickly run this script in high-mem interactive session:
# setwd("/home/groups/manishad/SAPH"); source("stitch_on_sherlock_SAPH.R")

# # load command line arguments
# args = commandArgs(trailingOnly = TRUE)
# start.num = as.numeric( args[1] )  # starting results number to stitch
# stop.num = as.numeric( args[2] )  # stopping results number to stitch



path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")

# PRELIMINARIES ----------------------------------------------

library(data.table)
library(dplyr)
library(testthat)
# s = stitch_files(.results.singles.path = "/home/groups/manishad/SAPH/sim_results/long",
#                  .results.stitched.write.path = "/home/groups/manishad/SAPH/sim_results/overall_stitched",
#                  .name.prefix = "long_results",
#                  .stitch.file.name="stitched.csv")

.results.singles.path = "/home/groups/manishad/SAPH/long_results"
.results.stitched.write.path = "/home/groups/manishad/SAPH/overall_stitched"
.name.prefix = "long_results"
.stitch.file.name="stitched.csv"


# MAKE STITCHED DATA ----------------------------------------------

# get list of all files in folder
all.files = list.files(.results.singles.path, full.names=TRUE)

# we only want the ones whose name includes .name.prefix
keepers = all.files[ grep( .name.prefix, all.files ) ]
length(keepers)

# grab variable names from first file
names = names( read.csv(keepers[1] ) )

# read in and rbind the keepers
tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )

# sanity check: do all files have the same names?
# if not, could be because some jobs were killed early so didn't get doParallelTime
#  variable added at the end
#  can be verified by looking at out-file for a job without name "doParallelTime"
allNames = lapply( tables, names )
# # find out which jobs had wrong number of names
# lapply( allNames, function(x) all.equal(x, names ) )table(s$method)
# allNames[[1]][ !allNames[[1]] %in% allNames[[111]] ]

# bind_rows works even if datasets have different names
#  will fill in NAs
s <- do.call(bind_rows, tables)

names(s) = names( read.csv(keepers[1], header= TRUE) )

if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
# write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )

cat("\n\n nrow(s) =", nrow(s))
cat("\n nuni(s$scen.name) =", nuni(s$scen.name) )




# ~ Check for Bad Column Names ---------------------------

# not sure why this is needed - has NA columns at end
names(s)
any(is.na(names(s)))

if ( any(is.na(names(s))) ) {
  NA.names = which( is.na(names(s) ) )
  s = s[ , -NA.names ]
  
}

s = s %>% filter(!is.na(scen.name))


# check runtimes - HOURS
summary(s$doParallel.seconds/60^2)

# 90th quantile of HOURS needed (1 rep)
# mathur WITHOUT robma: 0.78
quantile(s$doParallel.seconds/60^2, probs = 0.95)


# ~ Write stitched.csv ---------------------------

setwd(.results.stitched.write.path)
fwrite(s, .stitch.file.name)
# 
# # also make a zipped version
string = paste("zip -m stitched.zip", .stitch.file.name)
system(string)


# LOOK FOR MISSED JOBS ----------------------------------------------

# run in Sherlock ml load R
path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

missed.nums = sbatch_not_run( "/home/groups/manishad/SAPH/long_results",
                              "/home/groups/manishad/SAPH/long_results",
                              .name.prefix = "long_results",
                              .max.sbatch.num = 400 )


if ( FALSE ){
  
  
  setwd( paste(path, "/sbatch_files", sep="") )
  for (i in missed.nums) {
    system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
  }
  
}



# ~ Optional: Quick Summary ---------------------------


table(s$rep.name)

# reps per scen
# should be equal to reps.per.scen / reps.in.doParallel
s %>% group_by(scen.name, method) %>%
  summarise(n())

as.data.frame( s %>% group_by(method, scen.name, k.pub.nonaffirm) %>%
                  summarise(n()) )




#### General summary
# increase width of console print area for df viewing joy
options("width"=200)

t = s %>% group_by(hack, method, k.pub.nonaffirm) %>%
  #filter(rep.name == 1) %>% # TEMP - keep only first rep
  #filter(rep.name == 1) %>% # TEMP - keep only first rep
  filter(method == "rtma-pkg") %>%
  #filter(method == "2psm") %>%
  #filter(method %in% c("rtma-pkg", "jeffreys-mcmc-max-lp-iterate") ) %>%
  summarise( reps = n(),
             EstFail = mean(is.na(Mhat)),
             #Mhat = meanNA(Mhat),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             MLo = meanNA(MLo),
             MHi = meanNA(MHi),
             # Shat = meanNA(Shat),
             MhatNA = mean(is.na(Mhat))
             #MhatRhatGt1.05 = mean(MhatRhat>1.05),
             #MhatRhatGt1.02 = mean(MhatRhat>1.02)
  ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))
as.data.frame(t)


# Stefan: c.f. project log
t = s %>% group_by(hack, method, k.pub.nonaffirm) %>%
  filter(stringent.hack == TRUE &
           hack == "DV" & 
           k.pub.nonaffirm == 100 &
           strategy.stefan == "firstsig" & 
           alternative.stefan == "greater") %>%
#filter( method %in% c("rtma-pkg") ) %>%
  #filter(k.pub.nonaffirm ==10 & t2a == 0) %>%  # c.f. RSM_0
  #filter(rep.name > 1) %>% # TEMP - keep only >first rep
  #filter(rep.name == 1) %>% # TEMP - keep only first rep
  #filter(method == "rtma-pkg") %>%
  summarise( scen = scen.name[1],
             reps = n(),
             EstFail = mean(is.na(Mhat)),
             Mhat = meanNA(Mhat),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             MLo = meanNA(MLo),
             MHi = meanNA(MHi)
             # Shat = meanNA(Shat),
             # MhatNA = mean(is.na(Mhat)),
             # MhatRhatGt1.05 = mean(MhatRhat>1.05),
             # MhatRhatGt1.02 = mean(MhatRhat>1.02)
  ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))

as.data.frame(t)




t = s %>%
  #filter(scen.name == 1) %>%
  filter(hack == "affirm") %>%
  filter( method %in% c("jeffreys-mcmc-max-lp-iterate", "rtma-pkg", "robma") ) %>%
  #filter(k.pub.nonaffirm ==10 & t2a == 0) %>%  # c.f. RSM_0
  #filter(rep.name != 1) %>% # TEMP - keep only later reps
  #filter(rep.name == 1) %>% # TEMP - keep only first rep
  group_by(method) %>%
  summarise( reps = n(),
             EstFail = mean(is.na(Mhat)),
             Mhat = meanNA(Mhat),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             MLo = meanNA(MLo),
             MHi = meanNA(MHi)
             # Shat = meanNA(Shat),
             # MhatNA = mean(is.na(Mhat)),
             # MhatRhatGt1.05 = mean(MhatRhat>1.05),
             # MhatRhatGt1.02 = mean(MhatRhat>1.02)
  ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))

as.data.frame(t)
# when filtering to rep.name == 1, n() should be equal to reps.per.scen / reps.in.doParallel



#### For normal, large sims

# summarize scens that have run
param.vars.manip = c("hack", "k.pub.nonaffirm", "t2a")
library(tableone)
temp = s %>% filter(method == "naive")  # keep only 1 row per iterate
CreateCatTable(data = s, vars = param.vars.manip )


# increase width of console print area for df viewing joy
options("width"=200)

# General summary
t = s %>% group_by(hack, method) %>%
  #filter(k.pub.nonaffirm == 100 & hack == "DV") %>%
  #filter(k.pub.nonaffirm == 100 & hack == "subgroup") %>%  # choose sample size
  #filter(k.pub.nonaffirm == 30) %>%  # choose sample size
  #filter(method == "robma") %>%
  filter( method %in% c("robma", "jeffreys-mcmc-max-lp-iterate", "jeffreys-mcmc-pmed") ) %>%
  #filter(k.pub.nonaffirm >0 & t2a == 0) %>%
  summarise( reps = n(),
             EstFail = mean(is.na(Mhat)),
             Mhat = meanNA(Mhat),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             MLo = meanNA(MLo),
             MHi = meanNA(MHi)
             # Shat = meanNA(Shat),
             # MhatNA = mean(is.na(Mhat)),
             # MhatRhatGt1.05 = mean(MhatRhat>1.05),
             # MhatRhatGt1.02 = mean(MhatRhat>1.02)
  ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))

as.data.frame(t)


# compare to results in RSM_0 (mathur env)
#bm: this should look like the upper left of Fig 3 (close to 0 bias), but doesn't! Why??
t = s %>% group_by(hack, method) %>%
  #filter(k.pub.nonaffirm == 100 & hack == "DV") %>%
  #filter(k.pub.nonaffirm == 100 & hack == "subgroup") %>%  # choose sample size
  #filter(k.pub.nonaffirm == 30) %>%  # choose sample size
  #filter(method == "robma") %>%
  filter( method %in% c("rtma-pkg")) %>%
  #filter( method %in% c("naive", "2psm")) %>%
  filter(k.pub.nonaffirm == 10 &
           #hack == "favor-best-affirm-wch" & # favor first affirm
           hack == "affirm" &
           t2a == 0^2) %>%
  summarise( scens = length(unique(scen.name)),
             reps = n(),
             EstFail = mean(is.na(Mhat)),
             Mhat = meanNA(Mhat),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             MLo = meanNA(MLo),
             MHi = meanNA(MHi),
             # Shat = meanNA(Shat),
             # MhatNA = mean(is.na(Mhat)),
             MhatRhatGt1.05 = meanNA(MhatRhat>1.05),
             MhatRhatGt1.02 = meanNA(MhatRhat>1.02)
  ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))

as.data.frame(t)

# sanity check
temp = s[ s$k.pub.nonaffirm == 10 &
     s$hack == "affirm" & # favor first affirm
     s$t2a == 0^2 &
       s$method == "jeffreys-mcmc-max-lp-iterate", ]
meanNA(temp$Mhat)
median(temp$Mhat, na.rm = TRUE)
meanNA(temp$Mhat - temp$Mu)
meanNA(temp$MLo < temp$Mu & temp$MHi > temp$Mu)

meanNA(temp$MhatRhat > 1.05)

# before increasing stan parameters:
# as.data.frame(t)
#     hack                       method scens reps EstFail Mhat MhatBias MhatCover MhatWidth   MLo  MHi
# 1 affirm jeffreys-mcmc-max-lp-iterate     1  500       0 0.16    -0.34      0.98      4.76 -0.15 4.60
# 2 affirm           jeffreys-mcmc-pmed     1  500       0 0.31    -0.19      0.98      4.76 -0.15 4.60
# 3 affirm                        robma     1  500       0 0.30    -0.20      0.98      0.52  0.06 0.58




# summarize the scens that have run - stefan

t = s %>% group_by(hack, alternative.stefan, strategy.stefan) %>%
  summarise(n()) 
data.frame(t)
nrow(t)

table(s$hack)
table(s$alternative.stefan, s$strategy.stefan)
table(s$strategy.stefan)
table(s$k.pub.nonaffirm)

table(s$hack, s$k.pub.nonaffirm)





#** temp - explore robma issues
x = s %>% filter(k.pub.nonaffirm == 100 & hack == "DV") %>%  # choose sample size
  filter(method == "robma")

mean(x$MLo > x$Mhat)
mean(x$MHi < x$Mhat)

# by individual scen params
t = s %>% group_by(scen.name,
                   k.pub.nonaffirm,
                   hack,
                   strategy.stefan,
                   alternative.stefan,
                   method) %>%

  summarise( reps = n(),
             EstFail = mean(is.na(Mhat)),
             Mhat = meanNA(Mhat),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo)
             # MLo = meanNA(MLo),
             # MHi = meanNA(MHi),
             #Shat = meanNA(Shat),
             # MhatRhatGt1.05 = mean(MhatRhat>1.05),
             # MhatRhatGt1.02 = mean(MhatRhat>1.02)
             ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))

as.data.frame(t)



t = s %>% group_by(scen.name, k.pub.nonaffirm, Mu,
                   t2a, t2w, 
                   true.sei.expr,
                   method) %>%
  filter( method == "jeffreys-mcmc-max-lp-iterate" ) %>%
  # filter( grepl("jeffreys-mcmc", method) ) %>%
  summarise( reps = n(),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             # MLo = meanNA(MLo),
             # MHi = meanNA(MHi),
             Shat = meanNA(Shat),
             MhatNA = mean(is.na(Mhat)),
             MhatRhatGt1.05 = mean(MhatRhat>1.05),
  MhatRhatGt1.02 = mean(MhatRhat>1.02)) %>%
  #filter(reps > 1000) %>%
    mutate_if(is.numeric, function(x) round(x,2))

as.data.frame(t)


# iterates with small Rhat
t = s %>% group_by(scen.name, k.pub.nonaffirm, Mu,
                   t2a, t2w, 
                   true.sei.expr,
                   method) %>%
  filter( method == "jeffreys-mcmc-max-lp-iterate" ) %>% &
            MhatRhat < 1.02) %>%
  # filter( grepl("jeffreys-mcmc", method) ) %>%
  summarise( reps = n(),
             MhatBias = meanNA(Mhat - Mu),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             # MLo = meanNA(MLo),
             # MHi = meanNA(MHi),
             # Shat = meanNA(Shat),
             MhatNA = mean(is.na(Mhat)) ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))


as.data.frame(t)


# max-lp iterate 
t = s %>% group_by(scen.name, k.pub.nonaffirm, Mu, t2a, t2w, method) %>%
  filter( method == "jeffreys-mcmc-pmed" &
            MhatRhat < 1.01) %>%
  summarise( reps = n(),
             MhatMn = meanNA(Mhat),
             MhatCover = meanNA(MLo < Mu & MHi > Mu),
             MhatWidth = meanNA(MHi - MLo),
             # MLo = meanNA(MLo),
             # MHi = meanNA(MHi),
             # Shat = meanNA(Shat),
             MhatNA = mean(is.na(Mhat)) ) %>%
  #filter(reps > 1000) %>%
  mutate_if(is.numeric, function(x) round(x,2))


as.data.frame(t)


# MAKE AGG DATA ----------------------------------------------

path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")

# if this says "problem with column OptimConverged", 
#  you just need to comment out the optim columns in make_agg_data
#  because you didn't run those methods
agg = make_agg_data(s)

setwd(.results.stitched.write.path)
fwrite(agg, "agg.csv")

cat("\n\n nrow(agg) =", nrow(agg))
cat("\n nuni(agg$scen.name) =", nuni(agg$scen.name) )

# OPTIONAL: agg that keeps only iterates with very low Rhat
agg2 = make_agg_data( s %>% filter(MhatRhat<1.02 & ShatRhat<1.02))
fwrite(agg2, "agg_RhatsLt1.02.csv")





# look again at failures
t = agg %>% group_by(k.pub.nonaffirm, method) %>%
  summarise( mean(MhatEstFail),
             mean(MhatCIFail),
             mean(MhatTestReject)
             #meanNA(OptimxNAgreeOfConvergersMhatWinner)
             )
as.data.frame(t)


# errors of 2PSM when it fails
table( s$overall.error[ s$method == "2psm" & is.na(s$Mhat) ] )



##### Move to Local #####

# # stitched and agg -> local directory
# scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/overall_stitched/* /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Linked\ to\ OSF\ \(SAPH\)/Sherlock\ simulation\ results/Pilot\ simulations


