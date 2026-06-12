#!/bin/bash

#change time, mem, and cpu as appropriate to batch jobs submissions

#SBATCH --job-name=02_run_archr
#SBATCH --account=def-itobias
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt

#create positional parameter to input script
INPUT_RSCRIPT=${1}

#load module
module load apptainer/1.4.5

#run rscript
apptainer exec archr_image.sif Rscript $INPUT_RSCRIPT
