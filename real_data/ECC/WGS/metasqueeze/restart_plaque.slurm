#!/bin/bash
# The interpreter used to execute the script
# restart the squeezemeta pipeline since the initial run time out
#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=restart_plaque
#SBATCH --mail-type=BEGIN,END,ERROR
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128g
#SBATCH --time=48:00:00
#SBATCH --account=ligen0
#SBATCH --partition=standard
#SBATCH --output=plaque_metasqueeze_restart_%a.log
#SBATCH --error=plaque_metasqueeze_restart_error_%a.log


source /home/${USER}/.bashrc
micromamba activate SqueezeMeta
#filename=$(printf 'saliva_postHR_%02d.txt' ${SLURM_ARRAY_TASK_ID})

## output dir
WGS_folder=/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/WGS
cd ${WGS_folder}/metasqueeze


reads=${WGS_folder}/seqfetch/plaque_sequences
samplefile=${WGS_folder}/metasqueeze/plaque_list.txt
projectname=plaque_assembly

SqueezeMeta.pl --restart  -p $projectname  -step 10 --force_overwrite
