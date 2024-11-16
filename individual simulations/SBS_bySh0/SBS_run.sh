#!/bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=SBS_Multiple_Jobs
#SBATCH --output=output_%A_%a.out
#SBATCH --error=error_%A_%a.err
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --array=1-30                # Array range (submit 30 jobs)



# Execute your program, using SLURM_ARRAY_TASK_ID to differentiate between jobs
/flash/KondrashovU/FedyaK/SBS_bySh0/target/release/SBS_bySh0 --output output_${SLURM_ARRAY_TASK_ID}.txt


