#!/bin/sh

#SBATCH --job-name=ECC_WGS_locom
#SBATCH --time=04:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-4
#SBATCH --mem=90g
#SBATCH --cpus-per-task=6
#SBATCH --account=ligen0
#SBATCH --partition=standard
#SBATCH --export=ALL
#SBATCH --output=/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data/ECC/WGS/logs/%x/%x-%j.out
#SBATCH --error=/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data/ECC/WGS/logs/%x/%x-%j-error.out

module load gcc
module load fftw
module load gsl
module load Rtidyverse/4.2.0
Rscript --vanilla /scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data/ECC/WGS/ECCWGS_locom.R
