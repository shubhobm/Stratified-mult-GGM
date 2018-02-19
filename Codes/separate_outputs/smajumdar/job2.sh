#!/bin/bash
#SBATCH  --account=statistics
#SBATCH --qos=statistics
#SBATCH --job-name=R_test   #Job name	
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=2000mb
#SBATCH --mail-user=smajumdar@ufl.edu

#Record the time and compute node the job ran on
date; hostname; pwd

#Use modules to load the environment for R
module load R

#Run R script 
nohup nice -n19 R CMD BATCH sim_sep_1.R
