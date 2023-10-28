library(dacomp)
library(tidyverse)
library(splitstackshape)

rm(list=ls())
WGS_folder <- '/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/WGS/processed_data'
output_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data/ECC/WGS"
metadata_folder <- '/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/metadata'

load(file.path(WGS_folder, "Humann3.Rdata"))
load(file.path(WGS_folder, "kegg_and_taxa_slim.Rdata"))
load(file.path(metadata_folder, "metag_metashort.Rdata"))


fileid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

abund_mat <- switch(fileid,
                    plaque_kegg_abund,
                    saliva_kegg_abund,
                    plaque_species_abund,
                    saliva_species_abund)

dataname <- switch(fileid,
                   "Plaque_KEGG",
                   "Saliva_KEGG",
                   "Plaque_Species",
                   "Saliva_Species")


KEGGs_list=list('Plaque_KEGG'=plaque_kegg_paths, 'Saliva_KEGG'=saliva_kegg_paths)

prevalence <- rowMeans(abund_mat != 0)
gene_filter <- prevalence > 0.2

abund_mat_filter <- abund_mat[gene_filter, ]
#make new dataframe with samples in same order
metadata=metag_metashort[match(colnames(abund_mat_filter), metag_metashort$SampleName), ]
#change rownames to sample names
rownames(metadata)=metadata$SampleName
#check order 
if(all(rownames(metadata)==colnames(abund_mat_filter))!=T){stop('Rownames in wrong order')}

metadata$iscase <- metadata$CaseStatus == "case"

dacomp.ref.otus = dacomp.select_references(X = t(abund_mat_filter), median_SD_threshold = 1.3, verbose = F)
# test.method = ifelse(length(unique(Y)) == 2, DACOMP.TEST.NAME.WILCOXON, DACOMP.TEST.NAME.SPEARMAN)
res.dacomp <-  dacomp.test(X = t(abund_mat_filter), y = metadata$iscase, 
                          ind_reference_taxa = dacomp.ref.otus, 
                          test = DACOMP.TEST.NAME.WILCOXON, verbose = F, q = 0.05)
filename <- sprintf("%s_dacomp.rds", dataname)

saveRDS(res.dacomp, file.path(output_folder, filename))

