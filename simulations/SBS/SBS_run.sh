#!/bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=SBS_Multiple_Jobs
#SBATCH --output=output_%A_%a.out
#SBATCH --error=error_%A_%a.err
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# Path to the file containing the parameters
PARAM_FILE="parameters.txt"

# Read the specific line corresponding to the current array job index
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAM_FILE)

# Run the SBS program with the parameters read from the file
/flash/KondrashovU/FedyaK/SBS/target/release/SBS $PARAMS

