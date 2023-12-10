#!/bin/sh

#SBATCH --job-name=Dimension_experiment
#SBATCH --time=05:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-500
#SBATCH --mem=20g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out

module load gcc
module load fftw
module load gsl
module load Rtidyverse/4.2.0

cd $SLURM_SUBMIT_DIR

choice_array=({1..24..1})

for str in ${choice_array[@]}; do
  Rscript --vanilla simulate_Dimension.R -c ${str} -s ${SLURM_ARRAY_TASK_ID}
done


