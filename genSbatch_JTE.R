
# SET SIMULATION PARAMETERS MATRIX -----------------------------------------

# FOR CLUSTER USE
path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

allPackages = c("here",
                "magrittr",
                "dplyr",
                "data.table",
                # "fribidi",  # new dependency of tidyverse
                # "tidyverse", # these two can't be installed for some reason??
                "tidyr",
                "tibble",
                "metafor",
                "robumeta",
                "testthat",
                "truncdist",
                "gmm",
                "tmvtnorm",
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



# ### 2023-06-09 and 2023-06-11 - SIM.ENV = STEFAN ###
# scen.params = tidyr::expand_grid(
#   # without robma:
#   #rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; pet-peese ; jeffreys-mcmc ; rtma-pkg ; prereg-naive",
#   rep.methods = "robma",
#   
#   sim.env = "stefan",
# 
#   ### args shared between sim environments
#   #k.pub.nonaffirm = c(10, 15, 20, 30, 50, 70, 100),
#   k.pub.nonaffirm = c(10, 100, 30, 20),  # intentionally out of order so that jobs with boundary choices with complete first
#   hack = c("DV", "optstop", "subgroup"),
#   prob.hacked = c(0.8),
#   # important: if sim.env = stefan, these t2 args are ONLY used for setting start values
#   #   and for checking bias of Shat, so set them to have the correct t2a
#   #   not clear what t2w should be given the way stefan implements hacking
#   t2a = c(1),
#   t2w = c(0),
#   # same with Mu
#   Mu = c(0),
# 
#   # ### only needed if sim.env = "mathur": args from sim_meta_2
#   # Nmax = 30,
#   # m = 50,
#   #
#   # true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
#   # rho = c(0),
#   # ### end of stuff for sim.env = "mathur"
# 
#   ### only needed if sim.env = "stefan": args from sim_meta_2
#   strategy.stefan = c("firstsig", "smallest"),  # "firstsig" or "smallest"
#   alternative.stefan = c("greater", "two.sided"),  # "two.sided" or "greater"
#   stringent.hack = TRUE,  # mathur sims always effectively use stringent.hack = TRUE
#   ### end of stuff for sim.env = "stefan"
# 
#   # Stan control args
#   stan.maxtreedepth = 25,
#   stan.adapt_delta = 0.995,
# 
#   get.CIs = TRUE,
#   run.optimx = FALSE )
# 
# # hack.type = optstop must have strategy.stefan = "firstsig"
# scen.params = scen.params[ !(scen.params$hack == "optstop" & scen.params$strategy.stefan == "smallest"), ]
# 
# table(scen.params$hack, scen.params$strategy.stefan)



### 2023-05-31 and 2023-06-12 - SIM.ENV = MATHUR ###

scen.params = tidyr::expand_grid(
  # without robma:
  rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; pet-peese ; rtma-pkg ; prereg-naive",
  #rep.methods = "robma",

  sim.env = "mathur",

  # *If you reorder the args, need to adjust wrangle_agg_local
  ### args shared between sim environments
  #k.pub.nonaffirm = c(10, 100, 50, 20, 30, 15, 70),  # intentionally out of order so that jobs with boundary choices with complete first
  k.pub.nonaffirm = c(10, 100, 30, 20),
  hack = c("affirm", "favor-best-affirm-wch", "affirm2"),
  prob.hacked = c(0.8),
  # important: if sim.env = stefan, these t2 args are ONLY used for setting start values
  #   and for checking bias of Shat, so set them to have the correct t2a
  #   not clear what t2w should be given the way stefan implements hacking
  t2a = c(0, 0.2^2, 0.3^2, 0.5^2),
  t2w = c(0.2^2),
  # same with Mu
  Mu = c(0.5),

  ### only needed if sim.env = "mathur": args from sim_meta_2
  Nmax = 30,
  m = 50,
  true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
  rho = c(0),
  ### end of stuff for sim.env = "mathur"

  # ### only needed if sim.env = "stefan": args from sim_meta_2
  # strategy.stefan = c("firstsig", "smallest"),  # "firstsig" or "smallest"
  # alternative.stefan = c("greater", "two.sided"),  # "two.sided" or "greater"
  # stringent.hack = TRUE,  # mathur sims always effectively use stringent.hack = TRUE
  # ### end of stuff for sim.env = "stefan"

  # Stan control args
  stan.maxtreedepth = 25,
  stan.adapt_delta = 0.995,

  get.CIs = TRUE,
  run.optimx = FALSE )



# ### RSM_0 VERSION - AS IN 2022-5-17 SIMS ###
# scen.params = tidyr::expand_grid(
#   # full list (save):
#   # rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; jeffreys-mcmc ; jeffreys-sd ; jeffreys-var ; mle-sd ; mle-var ; csm-mle-sd ; 2psm-csm-dataset ; prereg-naive",
#    rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; jeffreys-mcmc ; jeffreys-sd ; prereg-naive",
#   #rep.methods = "naive ; jeffreys-mcmc ; jeffreys-sd",
# 
#   # args from sim_meta_2
#   Nmax = 30,
#   Mu = c(0.5),
#   t2a = c(0, 0.2^2, 0.3^2, 0.5^2),
#   t2w = c(0.2^2),
#   m = 50,
# 
#   hack = c("favor-best-affirm-wch", "affirm", "affirm2"),
#   rho = c(0),
#   k.pub.nonaffirm = c(10, 15, 20, 30, 50, 70, 100),
#   prob.hacked = c(0.8),
# 
#   true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
#    # true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)",
#    #                   "rbeta(n = 1, 2, 5)",
#    #                   "draw_lodder_se()"),
# 
#   # Stan control args
#   #@INCREASED 2022-4-26
#   stan.maxtreedepth = 25,
#   stan.adapt_delta = 0.995,
# 
#   get.CIs = TRUE,
#   run.optimx = FALSE )



# # OLD - Do I still need these?
# # hold constant the number of UNHACKED studies to 800
# scen.params$k[ scen.params$k.hacked == 800 ] = 800*2
# table(scen.params$k, scen.params$k.hacked)
# 
# # remove nonsense combinations
# # rho > 0 is pointless if there's only 1 draw
# scen.params = scen.params %>% dplyr::filter( !(rho > 0 & Nmax == 1) )


#@FOR CLARITY, maybe set to NA any scen params that aren't used based on stefan vs. mathur?
# also don't include all the extra combos

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
path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

scen.params = fread("scen_params.csv")
( n.scen = nrow(scen.params) )


# number of sbatches to generate (i.e., iterations within each scenario)
n.reps.per.scen = 500  
# ~ *** set sim.reps  -------------------------------------------------
n.reps.in.doParallel = 25
( n.files = ( n.reps.per.scen / n.reps.in.doParallel ) * n.scen )




scen.name = rep( scen.params$scen, each = ( n.files / n.scen ) )
jobname = paste("job", 1:n.files, sep="_")
# for re-run only: jobname = paste("job", 100:147, sep="_")
outfile = paste("/home/groups/manishad/SAPH/rmfiles/rm_", 1:n.files, ".out", sep="")
errorfile = paste("/home/groups/manishad/SAPH/rmfiles/rm_", 1:n.files, ".err", sep="")
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
                            jobtime = "02:00:00",
                            quality = "normal",
                            node_number = 1,
                            mem_per_node = 64000,
                            mailtype =  "NONE",
                            user_email = "mmathur@stanford.edu",
                            tasks_per_node = 16,
                            cpus_per_task = 1,
                            path_to_r_script = paste(path, "/doParallel_SAPH.R", sep=""),
                            args_to_r_script = paste("--args", jobname, scen.name, sep=" "),
                            write_path,
                            stringsAsFactors = F,
                            server_sbatch_path = NA)

generateSbatch(sbatch_params, runfile_path)

n.files

# run just the first one
#     sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/1.sbatch


# 2023-06-12 - 960 - mathur all other methods
# 2023-06-11 - 400 - stefan with only robma
# 2023-06-09 - 80 - stefan with all other methods
# 2023-05-30 - 480 - mathur with only robma
path = "/home/groups/manishad/SAPH"
setwd( paste(path, "/sbatch_files", sep="") )
for (i in 1:n.files) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
}



######## If Running Only Some Jobs To Fill Gaps ########

# run in Sherlock ml load R
path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

missed.nums = sbatch_not_run( "/home/groups/manishad/SAPH/long_results",
                              "/home/groups/manishad/SAPH/long_results",
                              .name.prefix = "long_results",
                              .max.sbatch.num = n.files )



setwd( paste(path, "/sbatch_files", sep="") )
for (i in missed.nums) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
}