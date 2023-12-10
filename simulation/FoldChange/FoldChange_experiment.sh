#!/bin/sh

#SBATCH --job-name=FoldChange_experiment
#SBATCH --time=03:00:00
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

choice_array=(1 2 3 4 5 6 7 8 9 10)

for str in ${choice_array[@]}; do
  Rscript --vanilla simulate_FoldChange.R -c ${str} -s ${SLURM_ARRAY_TASK_ID}
done

