
library(tidyverse)
library(DESeq2)
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
#convert to deseq 2 data
dds=DESeqDataSetFromMatrix(countData=abund_mat_filter, colData=metadata, design= ~ CaseStatus)

#refactor case status
dds$CaseStatus=factor(dds$CaseStatus, levels=c('control', 'case'))
#run deseq 2
dds2=DESeq(dds)
#pull results 
results=results(dds2, name='CaseStatus_case_vs_control')

if(fileid %in% c(1, 2)){
  #add kegg path annotation data
  KEGG_paths=as.data.frame(KEGGs_list[[dataname]])%>%rownames_to_column(var='KEGG')
  results_df=data.frame(results)%>%rownames_to_column(var='KEGG')%>%left_join(KEGG_paths)
  #split kegg annotations into columns
  colnames(results_df)=c(colnames(results_df)[1:7], 'PathAnnotation')
  r2=cSplit(results_df, 'PathAnnotation', '\\|', fixed=F, drop=F)
  r2=cSplit(r2, 'PathAnnotation_01', ';', fixed=F); r2=cSplit(r2, 'PathAnnotation_02', ';', fixed=F)
}
if(fileid %in% c(3, 4)){
  r2=data.frame(results)%>%rownames_to_column(var='Taxa')
}

filename <- sprintf("%s_deseq.rds", dataname)
saveRDS(r2, file.path(output_folder, filename))


