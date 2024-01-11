
# SET SIMULATION PARAMETERS MATRIX -----------------------------------------

# FOR CLUSTER USE
path = "/home/groups/manishad/JTE"
setwd(path)
source("helper_JTE.R")

allPackages = c("here",
                "magrittr",
                "dplyr",
                "data.table",
                # "fribidi",  # new dependency of tidyverse
                # "tidyverse", # these two can't be installed for some reason??
                "tidyr",
                "tibble",
                "testthat",
                "gmm",
                "doParallel",
                "foreach")
 

( packagesNeeded = allPackages[ !( allPackages %in% installed.packages()[,"Package"] ) ] )
if( length(packagesNeeded) > 0 ) install.packages(packagesNeeded)

# load all packages
lapply( allPackages,
        require,
        character.only = TRUE)

#**you need to see all "TRUE" printed by this in order for the package to actually be loaded

# set up sim params for cluster


# IMPORTANT NOTES ABOUT SCEN PARAMS:
# - Note that if you don't include any of these: jeffreys-sd ; jeffreys-var ; mle-sd ; mle-var
#  then you'll need to comment out Optim variables from the analysis.vars in make_agg_data and 
#  also from mutate in there
# - I think a similar thing will be true with the Rhats if you omit jeffreys-mcmc?
# - Usually good to run naive because it affects start values for subsequent methods (i.e., prevents
#   the start values from being the true ones)



### 2024-01-10 ###
scen.params = tidyr::expand_grid(
  # full list (save):
  rep.methods = "REML ; DL ; DL2 ; PM ; robu ; jeffreys",
  
  # *If you reorder the args, need to adjust wrangle_agg_local
  ### args shared between sim environments
  k.pub = c(100, 10,
            2, 3, 5, 20),  # intentionally out of order so that jobs with most interesting choices with complete first
 
  t2a = c(0.1^2,
          0.01^2, 0.05^2, 0.2^2, 0.5^2),

  # same with Mu
  Mu = c(0, 0.5, 1.1, 2.3), # for same as Langan's log-ORs
  true.dist = c("norm", "expo"),
  p0 = c(NA, 0.01, 0.05, 0.1, 0.5),  # as in Langan
  
  Ytype = c("bin-OR", "cont-SMD"),
  minN = c(40, 400, 2000),
  muN = c(40, 220, 400, 3000),
  
  # Stan control args
  stan.maxtreedepth = 25,
  stan.adapt_delta = 0.995,
  
  get.CIs = TRUE,
  run.optimx = FALSE )

table(scen.params$p0, useNA = "ifany")


#### Remove unwanted combinations

# ... of Mu and Ytype
remove = (scen.params$Mu != 0.5) & (scen.params$Ytype == "cont-SMD")
scen.params = scen.params[!remove,]
# sanity check:
table(scen.params$Mu, scen.params$Ytype)

# ... of Ytype and p0
remove = rep(FALSE, nrow(scen.params))
remove[ !is.na(scen.params$p0) & (scen.params$Ytype == "cont-SMD") ] = TRUE
remove[ is.na(scen.params$p0) & (scen.params$Ytype == "bin-OR") ] = TRUE
scen.params = scen.params[!remove,]
# sanity check:
table(scen.params$p0, scen.params$Ytype, useNA = "ifany")


# ... of minN and muN
# Lagnan does have one other version, which is not uniform
keep = rep(FALSE, nrow(scen.params))
keep[ scen.params$minN == 40 & scen.params$muN %in% c(40, 220) ] = TRUE
keep[ scen.params$minN == 400 & scen.params$muN %in% c(400) ] = TRUE
keep[ scen.params$minN == 2000 & scen.params$muN %in% c(3000) ] = TRUE
scen.params = scen.params[keep,]
table(scen.params$minN, scen.params$muN)





# add scen numbers
start.at = 1
scen.params = scen.params %>% add_column( scen = start.at : ( nrow(scen.params) + (start.at - 1) ),
                                          .before = 1 )


( n.scen = nrow(scen.params) )
# look at it
head( as.data.frame(scen.params) )

# write the csv file of params (to Sherlock)
setwd(path)
write.csv( scen.params, "scen_params.csv", row.names = FALSE )


########################### GENERATE SBATCHES ###########################

# load functions for generating sbatch files
path = "/home/groups/manishad/JTE"
setwd(path)
source("helper_JTE.R")

scen.params = fread("scen_params.csv")
( n.scen = nrow(scen.params) )


# number of sbatches to generate (i.e., iterations within each scenario)
n.reps.per.scen = 1000  
# ~ *** set sim.reps  -------------------------------------------------
n.reps.in.doParallel = 1000
( n.files = ( n.reps.per.scen / n.reps.in.doParallel ) * n.scen )




scen.name = rep( scen.params$scen, each = ( n.files / n.scen ) )
jobname = paste("job", 1:n.files, sep="_")
# for re-run only: jobname = paste("job", 100:147, sep="_")
outfile = paste("/home/groups/manishad/JTE/rmfiles/rm_", 1:n.files, ".out", sep="")
errorfile = paste("/home/groups/manishad/JTE/rmfiles/rm_", 1:n.files, ".err", sep="")
write_path = paste(path, "/sbatch_files/", 1:n.files, ".sbatch", sep="")
runfile_path = paste(path, "/testRunFile.R", sep="")


sbatch_params <- data.frame(jobname,
                            outfile,
                            errorfile,
                            # ma jobtimes by partition: sh_part
                            #jobtime = "02:00:00",  #@when running optimx methods, used sim.reps=100 and 5:00:00 here
                            
                            # 2023-05-29: mathur with only robma: first ran with 8:00, then with 24:00:00 to get the ~25% that timed out
                            # for RSM_1 sims with sim.env=stefan, n.reps.per.scen=500, and n.reps.in.doParallel=20 (1750 files):
                            # how to specify job times: https://www.sherlock.stanford.edu/docs/advanced-topics/job-management/#job-submission-limits
                            # days-hh:mm:ss
                            #jobtime = "1-00:00:00",  # 1 day
                            jobtime = "03:00:00",
                            quality = "normal",
                            node_number = 1,
                            mem_per_node = 64000,
                            mailtype =  "NONE",
                            user_email = "mmathur@stanford.edu",
                            tasks_per_node = 16,
                            cpus_per_task = 1,
                            path_to_r_script = paste(path, "/doParallel_JTE.R", sep=""),
                            args_to_r_script = paste("--args", jobname, scen.name, sep=" "),
                            write_path,
                            stringsAsFactors = F,
                            server_sbatch_path = NA)

generateSbatch(sbatch_params, runfile_path)

n.files

# run just the first one
#     sbatch -p qsu,owners,normal /home/groups/manishad/JTE/sbatch_files/1.sbatch


# 2024-01-10: 3120
path = "/home/groups/manishad/JTE"
setwd( paste(path, "/sbatch_files", sep="") )
for (i in 2:1000) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/JTE/sbatch_files/", i, ".sbatch", sep="") )
}



######## If Running Only Some Jobs To Fill Gaps ########

# run in Sherlock ml load R
path = "/home/groups/manishad/JTE"
setwd(path)
source("helper_JTE.R")

missed.nums = sbatch_not_run( "/home/groups/manishad/JTE/long_results",
                              "/home/groups/manishad/JTE/long_results",
                              .name.prefix = "long_results",
                              .max.sbatch.num = n.files )



setwd( paste(path, "/sbatch_files", sep="") )
for (i in missed.nums) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/JTE/sbatch_files/", i, ".sbatch", sep="") )
}