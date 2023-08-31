#!/bin/sh

#SBATCH --job-name=DirGamma_generation
#SBATCH --time=00:10:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=3g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=${SLURM_SUBMIT_DIR}/logs/%x-%a.out
#SBATCH --error=${SLURM_SUBMIT_DIR}/logs/%x-%a-error.out

module load gcc
module load fftw
module load gsl
module load Rtidyverse/4.2.0

cd $SLURM_SUBMIT_DIR
Rscript --vanilla global_gut_simulation.R -n 100 -p 0.1 -m mix -s ${SLURM_ARRAY_TASK_ID}
Rscript --vanilla global_gut_simulation.R -n 100 -p 0.1 -m enrich -s ${SLURM_ARRAY_TASK_ID}
Rscript --vanilla global_gut_simulation.R -n 100 -p 0.1 -m deplete -s ${SLURM_ARRAY_TASK_ID}
