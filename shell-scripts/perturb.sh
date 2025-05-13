#!/bin/bash

#SBATCH --job-name=PERTURB      # Set the job name
#SBATCH --nodes 1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --gpus a100:1
#SBATCH --mem 32gb
#SBATCH --time 04:00:00

#load modules
module load anaconda3/2023.09-0

#Go to working directory
# cd gem-diff/GEMDiff/GEMDiff

#Activate the created conda environment
source activate GEMDiff  #create before running script

#Perturb samples
###Make sure to change the model path which is found in the train log file.###
python scripts/perturb.py --config data/UTERN-UCECT.yaml --dir log-degs  --model_path log-degs/2025-05-08-21-16/model10000.pt --valid 
