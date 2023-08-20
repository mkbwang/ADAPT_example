library(POLDA)
library(dispmod)
library(survival)
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

res.polda <- polda(otu_table=abund_mat_filter,
                            metadata=metadata,
                            covar="iscase", ratio_model="loglogistic")

refgene <- res.polda$Reference_Taxa
refcounts <- colSums(abund_mat_filter[refgene, ])
polda_df <- res.polda$P_Value
polda_df$glm_effect <- 0
polda_df$glm_pval <- 0
polda_df$glm_padj <- 0
polda_df$cox_effect <- 0
polda_df$cox_pval <- 0
polda_df$cox_padj <- 0
for (j in 1:nrow(abund_mat_filter)){
  
  singlegene_count <- abund_mat_filter[polda_df$Taxon[j], ]
  singlegene_ratios <- rep(0, length(singlegene_count))
  exist <- singlegene_count > 0
  singlegene_ratios[exist] <- refcounts[exist] / singlegene_count[exist] 
  singlegene_ratios[!exist] <- refcounts[!exist]+1

  
  logit_model <- glm(cbind(singlegene_count, refcounts) ~ metadata$iscase, family="binomial")
  overdispersed_logit_model <- glm.binomial.disp(logit_model, verbose=F) |> summary()
  polda_df$glm_effect[j] <- overdispersed_logit_model$coefficients[2, 1]
  polda_df$glm_pval[j] <- overdispersed_logit_model$coefficients[2, 4]
  
  cox_model <- coxph(Surv(singlegene_ratios, exist) ~ metadata$iscase) |> summary()
  polda_df$cox_effect[j] <- cox_model$coefficients[1, 1]
  polda_df$cox_pval[j] <- cox_model$coefficients[1, 5]
  
}


polda_df$glm_padj <- p.adjust(polda_df$glm_pval, method="BH")
polda_df$cox_padj <- p.adjust(polda_df$cox_pval, method="BH")

res.polda$P_Value <- polda_df

filename <- sprintf("%s_polda.rds", dataname)

saveRDS(res.polda, file.path(output_folder, filename))

