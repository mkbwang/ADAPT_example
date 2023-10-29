# metasqueeze input requires a list of samples names

rm(list=ls())
folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/WGS/"

WGS_saliva <- read.csv(file.path(folder, 'seqfetch', 'runinfo_WGS_saliva.csv'))

Runid_saliva <- rep(WGS_saliva$Run, each=2)
sample_saliva <- rep(WGS_saliva$Library.Name, each=2)
directions <- rep(c(1,2), 30)
saliva_filenames <- sprintf("%s_%d.fastq.gz", Runid_saliva, directions)
pair_names <- sprintf("pair%d", directions)

squeezelist_saliva <- data.frame(sample_saliva, saliva_filenames, pair_names)
write.table(squeezelist_saliva, file=file.path(folder, 'metasqueeze', 'saliva_list.txt'),
            sep='\t', row.names=F,col.names=F, quote=F)



WGS_plaque <- read.csv(file.path(folder, 'seqfetch', 'runinfo_WGS_plaque.csv'))

Runid_plaque <- rep(WGS_plaque$Run, each=2)
sample_plaque <- rep(WGS_plaque$Library.Name, each=2)
directions <- rep(c(1,2), 30)
plaque_filenames <- sprintf("%s_%d.fastq.gz", Runid_plaque, directions)
pair_names <- sprintf("pair%d", directions)

squeezelist_plaque <- data.frame(sample_plaque, plaque_filenames, pair_names)
write.table(squeezelist_plaque, file=file.path(folder, 'metasqueeze', 'plaque_list.txt'),
            sep='\t', row.names=F,col.names=F, quote=F)

