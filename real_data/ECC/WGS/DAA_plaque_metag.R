library(dplyr)
rm(list=ls())
folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/WGS"

phy_metag_plaque <- readRDS(file.path(folder, "phyloseq", "phy_metag_plaque.rds"))


count_table <- otu_table(phy_metag_plaque)@.Data
prevalences <- rowMeans(count_table != 0)
count_table_subset <- count_table[prevalences > 0.05, ]
taxa_names <- rownames(count_table_subset)
metadata$iscase <- metadata$CaseStatus == "case"

set.seed(1)
# ADAPT
library(ADAPT)

adapt_output_noboot <- adapt(otu_table=count_table_subset, metadata=metadata, covar="iscase",
                      adjust=NULL, prevalence_cutoff=0.05, genes_are_rows=T, boot=F)
adapt_output_boot <- adapt(otu_table=count_table_subset, metadata=metadata, covar="iscase",
                           adjust=NULL, prevalence_cutoff=0.05, genes_are_rows=T, boot=T, n_boot_gene=100,
                           boot_replicate=1000, alpha=0.05)
adapt_result_df <- adapt_output_boot$P_Value

adapt_reference_df <- adapt_result_df %>% filter(Gene %in% adapt_output_boot$Reference_Gene)
adapt_DA_df <- adapt_result_df %>% filter(Gene %in% adapt_output_boot$DA_Gene)

adapt_DAtaxa <- adapt_output_boot$DA_Gene


# Aldex
library(ALDEx2)
aldex_result_df <- aldex(reads=count_table_subset,
                         conditions=metadata$iscase)
aldex_DAtaxa <- rownames(aldex_result_df)[aldex_result_df$wi.eBH < 0.05]



# MetagenomeSeq
library(metagenomeSeq)
MRmetadata <- AnnotatedDataFrame(metadata)
MRobj <- newMRexperiment(counts=count_table_subset,
                         phenoData=MRmetadata)

MRobj <- wrenchNorm(MRobj, condition=metadata$CaseStatus)
mod <- model.matrix(~iscase, data=pData(MRobj))
MRfit <- fitFeatureModel(MRobj, mod)
MR_coefficients <- MRcoefs(MRfit, number=nrow(count_table_subset), adjustMethod="BH")
MR_DAtaxa <- rownames(MR_coefficients)[MR_coefficients$adjPvalues < 0.05] |> na.omit()


# DACOMP
library(dacomp)
set.seed(1)
selected_references <- dacomp.select_references(
  X = count_table_subset, verbose = F, run_in_parallel = T, Nr.Cores=4,
  minimal_TA=10)
updated_selected_references <- dacomp.validate_references(X = t(count_table_subset), #Counts
                                                          Y = metadata$iscase, #Traits
                                                          ref_obj = selected_references, #reference checked, must be object from dacomp.select_references(...)
                                                          test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                          Q_validation = 0.1, #FDR level for checking
                                                          Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                          Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                          NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                          Verbose = F)
dacomp_output <- dacomp.test(X = t(count_table_subset), #counts data
                             y = metadata$iscase, #phenotype in y argument
                             # obtained from dacomp.select_references(...):
                             ind_reference_taxa = updated_selected_references,
                             test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                             disable_DSFDR = T,
                             q=0.05,
                             verbose = F)

dacomp_DAtaxa <- taxa_names[dacomp_output$p.values.test.adjusted < 0.05] |>
  na.omit()


# RDB
library(RDB)

depths <- rowSums(t(count_table_subset))
compositions <- t(count_table_subset) / depths
rdb_output <- rdb(P=compositions, Z=metadata$iscase, alpha=0.05, fdr=T)

rdb_DAtaxa <- taxa_names[rdb_output]




# LOCOM
library(LOCOM)
locom_output <- locom(otu.table=t(count_table_subset),
                      Y=metadata$iscase,fdr.nominal=0.05, prev.cut=0.05,
                      n.rej.stop=100, seed=1, n.perm.max=2e4, Firth.thresh = 1,
                      n.cores=4)

locom_qval <- locom_output$q.otu
locom_DAtaxa <- locom_output$detected.otu



# ANCOMBC
library(ANCOMBC)
ancombc_output <- ancombc2(phy_metag_plaque, fix_formula='CaseStatus',
                           p_adj_method='BH',
                           group='CaseStatus', struc_zero=FALSE, prv_cut=0.05, n_cl=4)

ancombc_DAtaxa <- ancombc_output$res$taxon[ancombc_output$res$q_CaseEverControl < 0.05]



# LinDA
library(LinDA)
linda_output <- linda(otu.tab=count_table_subset,
                      meta=data.frame(metadata),
                      formula="~iscase",
                      n.cores=4)
linda_df <- linda_output$output$iscase
linda_DAtaxa <- rownames(linda_df)[linda_df$reject]


# aggregate outputs
allDAtaxa <- Reduce(union, list(adapt_DAtaxa, aldex_DAtaxa, ancombc_DAtaxa,
                                dacomp_DAtaxa, linda_DAtaxa, locom_DAtaxa,
                                rdb_DAtaxa, MR_DAtaxa))


output_df <- data.frame(Taxon = allDAtaxa)
output_df$ADAPT <- output_df$Taxon %in% adapt_DAtaxa
output_df$ALDEx2 <- output_df$Taxon %in% aldex_DAtaxa
output_df$DACOMP <- output_df$Taxon %in% dacomp_DAtaxa
output_df$ANCOMBC <- output_df$Taxon %in% ancombc_DAtaxa
output_df$metagenomeSeq <- output_df$Taxon %in% MR_DAtaxa
output_df$LinDA <- output_df$Taxon %in% linda_DAtaxa
output_df$LOCOM <- output_df$Taxon %in% locom_DAtaxa
output_df$RDB <- output_df$Taxon %in% rdb_DAtaxa

taxonomy_table <- tax_table(phy_metag_plaque)@.Data
selected_taxonomies <- taxonomy_table[output_df$Taxon, ] |> as.data.frame()
output_df <- cbind(output_df, selected_taxonomies)


write.csv(output_df, file.path(folder, "method_comparison_WGS_plaque.csv"),
          row.names=T)
write.csv(adapt_result_df, file.path(folder, "ADAPT_WGS_plaque_result.csv"),
          row.names=F)

