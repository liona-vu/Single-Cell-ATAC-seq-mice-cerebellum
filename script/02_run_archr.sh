#!/bin/bash

#SBATCH --job-name=run_archr
#SBATCH --account=def-itobias
#SBATCH --time=4:00:00
#SBATCH --mem=16384M
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1          #threading parellelization
##SBATCH --error=%x_%j_error.txt
##SBATCH --output=%x_%j_output.txt

#create positional parameters to input sif image and script
ARCHR_IMAGE=${1}
INPUT_RSCRIPT=${2}

#load module
module load apptainer/1.4.5
#module load r/4.5.0

#run rscript
apptainer exec $ARCHR_IMAGE Rscript $INPUT_RSCRIPT


#for an interactive session on the HPC, run this instead
#apptainer shell archr_image.sif

#then run 
#R
