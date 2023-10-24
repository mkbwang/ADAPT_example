#!/bin/sh

#SBATCH --job-name=sim_metagenomeseq
#SBATCH --time=00:50:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=3g
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/wangmk/MDAWG/POLDA_example/simulation/MetagenomeSeq/slurm-sim-metagenomeseq.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA_example/simulation/MetagenomeSeq/metagenomeseq_simulation.R

