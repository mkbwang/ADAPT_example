#!/bin/sh

#SBATCH --job-name=sim_POLDA_nobiascorrect
#SBATCH --time=03:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=1g
#SBATCH --cpus-per-task=3
#SBATCH --output=/home/wangmk/MDAWG/POLDA_example/simulation/POLDA_nobiascorrect/slurm-sim-POLDA.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA_example/simulation/POLDA_nobiascorrect/polda_simulation.R FALSE
