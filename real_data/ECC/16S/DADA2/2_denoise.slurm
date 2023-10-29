#!/bin/sh
#SBATCH --job-name=ECC1_denoise
#SBATCH --account=ligen0
#SBATCH --partition=standard
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3
#SBATCH --mem=20GB
#SBATCH --time=04:00:00
#SBATCH --output=denoise_output.log
#SBATCH --error=denoise_error.log


#load R

module load Rtidyverse/4.2.0

Rscript 2_denoise_maketable.R
