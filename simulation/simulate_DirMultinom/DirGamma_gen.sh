#!/bin/sh

#SBATCH --job-name=DirGamma_generation
#SBATCH --time=02:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=4g
#SBATCH --cpus-per-task=3
#SBATCH --output=/home/wangmk/MDAWG/POLDA_example/simulation/simulate_DirMultinom/slurm-sim.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA_example/simulation/simulate_DirMultinom/simulate_DirMultinom.R
