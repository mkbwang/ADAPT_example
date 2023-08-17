#!/bin/sh

#SBATCH --job-name=FPR_cox
#SBATCH --time=01:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=4g
#SBATCH --cpus-per-task=3
#SBATCH --account=ligen0
#SBATCH --partition=standard
#SBATCH --export=ALL
#SBATCH --output=/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation/FPR/permutation/logs/%x/%x-%j.out
#SBATCH --error=/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation/FPR/permutation/logs/%x/%x-%j-error.out


module load fftw
module load gsl
module load Rtidyverse/4.2.0
Rscript --vanilla /scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation/FPR/permutation/test_cox.R
