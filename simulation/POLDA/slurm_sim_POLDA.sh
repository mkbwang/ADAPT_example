#!/bin/sh

#SBATCH --job-name=sim_POLDA_biascorrect
#SBATCH --time=01:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=2g
#SBATCH --cpus-per-task=2
#SBATCH --output=/home/wangmk/MDAWG/POLDA_example/simulation/POLDA/slurm-sim-POLDA.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA_example/simulation/POLDA/polda_simulation.R TRUE
