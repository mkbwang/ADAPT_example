#!/bin/sh

#SBATCH --job-name=sim_PTDA
#SBATCH --time=02:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=8g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out

module load gcc
module load fftw
module load gsl
module load Rtidyverse/4.2.0

cd $SLURM_SUBMIT_DIR


## baseline
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000


## Different fold change
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m mix -f 2.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m enrich -f 2.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m mix -f 3.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m enrich -f 3.5 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m mix -f 4 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m enrich -f 4 -s ${SLURM_ARRAY_TASK_ID} -d 100000


## Different proportion of DA
Rscript --vanilla ptda_simulation.R -n 100 -p 0.05 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.05 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.2 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.2 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000


## Different sample size
Rscript --vanilla ptda_simulation.R -n 50 -p 0.1 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000
Rscript --vanilla ptda_simulation.R -n 50 -p 0.1 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 100000


## smaller sequencing depths
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m mix -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 20000
Rscript --vanilla ptda_simulation.R -n 100 -p 0.1 -m enrich -f 3 -s ${SLURM_ARRAY_TASK_ID} -d 20000



# #Rscript --vanilla polda_simulation.R -n 50 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a lognormal
# #Rscript --vanilla polda_simulation.R -n 50 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a lognormal
# #Rscript --vanilla polda_simulation.R -n 50 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a lognormal
# Rscript --vanilla polda_simulation.R -n 50 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000  -a lognormal
# Rscript --vanilla polda_simulation.R -n 50 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000  -a lognormal
# Rscript --vanilla polda_simulation.R -n 50 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -a lognormal
# Rscript --vanilla polda_simulation.R -n 50 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a lognormal
# Rscript --vanilla polda_simulation.R -n 50 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a lognormal
# Rscript --vanilla polda_simulation.R -n 50 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a lognormal
# 
# 
# 
# #Rscript --vanilla polda_simulation.R -n 100 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a lognormal
# #Rscript --vanilla polda_simulation.R -n 100 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a lognormal
# #Rscript --vanilla polda_simulation.R -n 100 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a lognormal
# Rscript --vanilla polda_simulation.R -n 100 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0.5 -a lognormal
# Rscript --vanilla polda_simulation.R -n 100 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0.5 -a lognormal
# Rscript --vanilla polda_simulation.R -n 100 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0.5 -a lognormal
# Rscript --vanilla polda_simulation.R -n 100 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a lognormal
# Rscript --vanilla polda_simulation.R -n 100 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a lognormal
# Rscript --vanilla polda_simulation.R -n 100 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a lognormal
# # 
# # 
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0.5 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0.5 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 0.5 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m mix -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m deplete -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a loglogistic
# # Rscript --vanilla polda_simulation.R -n 100 -p 0 -m enrich -z 1 -s ${SLURM_ARRAY_TASK_ID} -d 100000 -t 1 -a loglogistic
# # 
# # 
