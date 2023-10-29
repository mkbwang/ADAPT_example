#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --job-name=ECCdownload_plaqueWGS
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mem-per-cpu=1g
#SBATCH --account=ligen0
#SBATCH --partition=standard
#SBATCH --output=download_plaque.log
#SBATCH --error=download_plaque_error.log


module load Bioinformatics
module load gcc
module load sratoolkit

# sed 1d runinfo_Amplicon16S.csv | cut -d "," -f 1 > SRR.numbers.amplicon
cat SRR.numbers.plaque.WGS | parallel prefetch
cat SRR.numbers.plaque.WGS | parallel fastq-dump --gzip --split-3 --dumpbase --skip-technical -O plaque_sequences
