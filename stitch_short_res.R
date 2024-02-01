
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




# PRELIMINARIES ----------------------------------------------

library(data.table)
library(dplyr)
library(testthat)
library(doParallel)

num_cores = 16
registerDoParallel(cores = 16)

path = "/home/groups/manishad/JTE"
setwd(path)
source("helper_JTE.R")
source("analyze_sims_helper_JTE.R")


# MAKE AGG DATA - FROM PRE-AGGREGATED SHORT RESULTS ----------------------------------------------

# can be run in interactive session - :)!
.results.singles.path = "/home/groups/manishad/JTE/short_results"
.results.stitched.write.path = "/home/groups/manishad/JTE/overall_stitched"
.name.prefix = "short_results"
.stitch.file.name="aggo.csv"


# get list of all files in folder
all.files = list.files(.results.singles.path, full.names=TRUE)

# we only want the ones whose name includes .name.prefix
keepers = all.files[ grep( .name.prefix, all.files ) ]
length(keepers)

# grab variable names from first file
names = names( read.csv(keepers[1] ) )


# read in and rbind the keepers in parallel
rbind_tables <- function(files) {
  tables <- lapply(files, fread)
  rbindlist(tables, fill = TRUE)
}

split_files <- split(keepers, 1:length(keepers) %% num_cores)

tables <- foreach(i = 1:length(split_files), .packages = c("data.table")) %dopar% {
  rbind_tables(split_files[[i]])
}



# combine the results
agg <- rbindlist(tables, fill = TRUE)
length(unique(agg$scen.name))


# write it
setwd(.results.stitched.write.path)
fwrite(agg, .stitch.file.name)

