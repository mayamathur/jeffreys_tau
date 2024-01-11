#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=job_stitch
#################  
#a file for job output, you can check job progress
#SBATCH --output=/home/groups/manishad/JTE/rmfiles/rm_stitch.out
#################
# a file for errors from the job
#SBATCH --error=/home/groups/manishad/JTE/rmfiles/rm_stitch.err
#################
#time you think you need; default is one hour
#SBATCH --time=00:30:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#submit to both owners and normal partition
#SBATCH -p normal,owners,qsu
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB
#SBATCH --mem=64000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=NONE
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=mmathur@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=16
#################
#SBATCH --cpus-per-task=1
#now run normal batch commands

ml load v8
ml load R/4.1.2
R -f /home/groups/manishad/JTE/stitch_on_sherlock_JTE.R
