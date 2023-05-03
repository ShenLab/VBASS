#!/bin/bash 
#
#SBATCH --job-name=ASD
# Default in slurm
#SBATCH --mail-user=gz2294@cumc.columbia.edu
#SBATCH --mail-type=ALL
#SBATCH -t 72:0:0 # Request 5 hours run time
# Define threads and memeory
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10gb
#SBATCH --output=ASD.log
date;hostname;pwd

source /home/gz2294/.bashrc
/share/vault/Users/gz2294/miniconda3/envs/r4-base/bin/Rscript run.ASD.GO.R ${SLURM_ARRAY_TASK_ID}

date
