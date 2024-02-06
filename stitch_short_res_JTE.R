
# can be run in interactive session - :)!

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


# write it
setwd(.results.stitched.write.path)
fwrite(agg, .stitch.file.name)

length(unique(agg$scen.name))

