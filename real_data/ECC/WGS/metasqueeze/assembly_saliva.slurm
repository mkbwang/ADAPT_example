#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=assembly_saliva
#SBATCH --mail-type=BEGIN,END,ERROR
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128g
#SBATCH --time=48:00:00
#SBATCH --account=ligen0
#SBATCH --partition=standard
#SBATCH --output=saliva_metasqueeze_%a.log
#SBATCH --error=saliva_metasqueeze_error_%a.log


source /home/${USER}/.bashrc
micromamba activate SqueezeMeta

## output dir
WGS_folder=/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/WGS
cd ${WGS_folder}/metasqueeze

export http_proxy=http://proxy.arc-ts.umich.edu:3128/
export FTP_PROXY=http://proxy.arc-ts.umich.edu:3128/
export ftp_proxy=http://proxy.arc-ts.umich.edu:3128/
export HTTPS_PROXY=http://proxy.arc-ts.umich.edu:3128/
export https_proxy=http://proxy.arc-ts.umich.edu:3128/
export HTTP_PROXY=http://proxy.arc-ts.umich.edu:3128/

reads=${WGS_folder}/seqfetch/saliva_sequences
samplefile=${WGS_folder}/metasqueeze/saliva_list.txt
projectname=saliva_assembly

SqueezeMeta.pl -m coassembly -s $samplefile -p $projectname  -f $reads -c 1000 -b 16 --euk

