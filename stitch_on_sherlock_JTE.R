


# run this interactively in ml load R or via:
#   sbatch -p qsu,owners,normal /home/groups/manishad/JTE/job_stitch.sbatch
# sacct --name=job_stitch
# look at its out file:
# less /home/groups/manishad/JTE/rmfiles/rm_stitch.out
# less /home/groups/manishad/JTE/rmfiles/rm_stitch.err
# shift-g to jump to bottom of files

# for non-huge simulations, can often run this script interactively in a higher-memory
#  Sherlock session:
# ml load R/4.1.2
# srun --mem=32G --time=3:00:00 --pty bash
# R


# to be run by stitch.sbatch or manually
# To quickly run this script in high-mem interactive session:
# setwd("/home/groups/manishad/JTE"); source("stitch_on_sherlock_JTE.R")

# # load command line arguments
# args = commandArgs(trailingOnly = TRUE)
# start.num = as.numeric( args[1] )  # starting results number to stitch
# stop.num = as.numeric( args[2] )  # stopping results number to stitch


path = "/home/groups/manishad/JTE"
setwd(path)
source("helper_JTE.R")
source("analyze_sims_helper_JTE.R")

# PRELIMINARIES ----------------------------------------------

library(data.table)
library(dplyr)
library(testthat)
# s = stitch_files(.results.singles.path = "/home/groups/manishad/JTE/sim_results/long",
#                  .results.stitched.write.path = "/home/groups/manishad/JTE/sim_results/overall_stitched",
#                  .name.prefix = "long_results",
#                  .stitch.file.name="stitched.csv")


.results.singles.path = "/home/groups/manishad/JTE/long_results"
.results.stitched.write.path = "/home/groups/manishad/JTE/overall_stitched"
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
tables <- lapply( keepers, function(x) {
  
  y = tryCatch( read.csv(x, header = TRUE), error = function(e) NULL )
  
  y[[ "N.expr" ]] = as.character(y[[ "N.expr" ]] )  # only needed if it is just a number, because it turns into a double for certain datasets and then can't be concatenated with character ones
  y
  } )


cat("\n\nFinished reading in tables")

# sanity check: do all files have the same names?
# if not, could be because some jobs were killed early so didn't get doParallelTime
#  variable added at the end
allNames = lapply( tables, names )
# # find out which jobs had wrong number of names
# lapply( allNames, function(x) all.equal(x, names ) )table(s$method)
# allNames[[1]][ !allNames[[1]] %in% allNames[[111]] ]

# bind_rows works even if datasets have different names
#  will fill in NAs
s <- do.call(bind_rows, tables)

cat("\n\nFinished s <- do.call")

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

if (FALSE) {
  
  path = "/home/groups/manishad/JTE"
  setwd(path)
  source("helper_JTE.R")
  
  # this will write missed_nums.csv to the location specified by second arg
  missed.nums = sbatch_not_run( "/home/groups/manishad/JTE/long_results",
                                "/home/groups/manishad/JTE/overall_stitched",
                                .name.prefix = "long_results",
                                .max.sbatch.num = 2496 )
  
  
  
  setwd( paste(path, "/sbatch_files", sep="") )
  for (i in missed.nums) {
    system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/JTE/sbatch_files/", i, ".sbatch", sep="") )
  }
  
}




# ~ Optional: Quick Summary ---------------------------

if (FALSE) {

  # reps per scen
  # should be equal to reps.per.scen / reps.in.doParallel
  s %>% group_by(scen.name, method) %>%
    summarise(n())
  
  as.data.frame( s %>% group_by(method, scen.name, k.pub) %>%
                   summarise(n()) )
  
  
  # summarize sim params that have run so far
  library(tableone)
  vars = c("k.pub", "t2a", "Mu", "true.dist", "p0", "Ytype", "N.expr")
  CreateTableOne( dat = s,
                  vars = vars,
                  factorVars = vars )
  
  
  #### Look at sanity checks for binary Y
  # increase width of console print area for df viewing joy
  options("width"=200)
  
  t = s %>%
    filter(method == "REML") %>%
    filter(Ytype == "bin-OR") %>%
    group_by(scen.name, Mu, p0) %>%
    
    summarise( reps = n(),
               
               MhatBias = meanNA(Mhat - Mu),
              ShatBias = meanNA(Mhat - Mu),
               sancheck_mean_pY0 = meanNA(sancheck_mean_pY0),
               sancheck_mean_pY = meanNA(sancheck_mean_pY)
               
    ) %>%
    #filter(reps > 1000) %>%
    mutate_if(is.numeric, function(x) round(x,2))
  as.data.frame(t)
  
  
  #sum(s$k.pub == 10 & s$t2a == 0.2^2 & s$Mu == 0)
  
  
  # main results
  t = s %>% group_by(method) %>%
    #filter(Ytype == "cont-SMD") %>% 
    filter(Ytype == "bin-OR") %>% 

    summarise( reps = n(),
               # EstFail = mean(is.na(Mhat)),
               # #Mhat = meanNA(Mhat),
               
               MhatBias = meanNA(Mhat - Mu),
               MhatCover = meanNA(MLo <= Mu & MHi >= Mu),
               MhatWidth = meanNA(MHi - MLo),
               MhatMSE = meanNA( ( Mhat - Mu )^2 ),
               
               ShatBias = meanNA(Shat - sqrt(t2a)),
               ShatCover = meanNA(SLo <= sqrt(t2a) & SHi >= sqrt(t2a)),
               ShatWidth = meanNA(SHi - SLo),
               ShatMSE = meanNA( ( Shat - sqrt(t2a) )^2 ),
               # SLo = meanNA(SLo),
               # SHi = meanNA(SHi),
               # Shat = meanNA(Shat),
               #ShatNA = mean(is.na(Shat)),
               

               # MLo = meanNA(MLo),
               # MHi = meanNA(MHi),
               # Shat = meanNA(Shat),
               #MhatNA = mean(is.na(Mhat))
               #MhatRhatGt1.05 = mean(MhatRhat>1.05),
               #MhatRhatGt1.02 = mean(MhatRhat>1.02)
    ) %>%
    #filter(reps > 1000) %>%
    mutate_if(is.numeric, function(x) round(x,2))
  as.data.frame(t)
  
}




# MAKE AGG DATA ----------------------------------------------

path = "/home/groups/manishad/JTE"
setwd(path)
source("helper_JTE.R")
source("analyze_sims_helper_JTE.R")

# if this says "problem with column OptimConverged", 
#  you just need to comment out the optim columns in make_agg_data
#  because you didn't run those methods
agg = make_agg_data(s)

setwd(.results.stitched.write.path)
fwrite(agg, "aggo.csv")

cat("\n\n nrow(agg) =", nrow(agg))
cat("\n nuni(agg$scen.name) =", nuni(agg$scen.name) )



##### Move to Local #####

# # stitched and agg -> local directory
# scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/JTE/overall_stitched/* /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(JTE\)/Linked\ to\ OSF\ \(JTE\)/Sherlock\ simulation\ results/Pilot\ simulations


