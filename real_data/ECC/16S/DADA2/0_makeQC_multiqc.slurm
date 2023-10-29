#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=make_multi_qc_ECCPaper1
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=BEGIN,END,ERROR
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB 
#SBATCH --time=10:00:00
#SBATCH --account=ligen0
#SBATCH --partition=standard
#SBATCH --output=qc_output.log
#SBATCH --error=qc_error.log

source /home/wangmk/.bashrc
micromamba activate multiqc # Mukai's local environment
module load Bioinformatics
module load fastqc

folder=/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/16S/seqfetch

## forward
cd ${folder}/raw_sequences/FWD
fastqc *.fastq.gz --quiet --outdir=${folder}/raw_sequences/FWD/fastqc/
multiqc ./ --quiet --force --filename  FWD_16S_multiqc_report.html  --outdir=${folder}/raw_sequences/FWD/fastqc/


## backward
cd ${folder}/raw_sequences/REV
fastqc *.fastq.gz --quiet  --outdir=${folder}/raw_sequences/REV/fastqc/
multiqc ./ --quiet --force --filename  REV_16S_multiqc_report.html --outdir=${folder}/raw_sequences/REV/fastqc/



