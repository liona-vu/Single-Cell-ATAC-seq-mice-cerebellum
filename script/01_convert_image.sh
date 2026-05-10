#!/bin/bash

#SBATCH --job-name=download_archr
#SBATCH --account=def-itobias
#SBATCH --time=00:10:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt

#load module
module load apptainer/1.4.5

#convert archr docker image to sif image for apptainer
apptainer pull docker://immanuelazn/archr:latest
