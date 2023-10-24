#!/bin/sh

#SBATCH --job-name=sim_LOCOM
#SBATCH --time=4:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --account=ligen0
#SBATCH --mem=8g
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out

module load gcc
module load fftw
module load gsl
module load Rtidyverse/4.2.0

## baseline
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000


## Different fold change
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m mix -f 2.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m enrich -f 2.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m mix -f 3.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m enrich -f 3.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m mix -f 4 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m enrich -f 4 -s ${SLURM_ARRAY_TASK_ID} -d 100000



## Different proportion of DA
Rscript --vanilla locom_simulations.R -n 100 -p 0.05 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.05 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.2 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 100 -p 0.2 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000


## Different sample size
Rscript --vanilla locom_simulations.R -n 50 -p 0.1 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla locom_simulations.R -n 50 -p 0.1 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000


## smaller sequencing depths
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 20000
Rscript --vanilla locom_simulations.R -n 100 -p 0.1 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 20000

