library(dplyr)
library(stringr)
rm(list=ls())
folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/WGS"

phy_metag_plaque <- readRDS(file.path(folder, "phyloseq", "phy_metag_plaque.rds"))

phy_metag_plaque <- phy_metag_plaque %>% filter_taxa(function(ct) mean(ct > 0) > 0.05, TRUE)
count_table <- otu_table(phy_metag_plaque)@.Data
# prevalences <- rowMeans(count_table != 0)
# count_table_subset <- count_table[prevalences > 0.05, ]
taxa_names <- rownames(count_table)
metadata <- data.frame(sample_data(phy_metag_plaque))
metadata$iscase <- metadata$CaseStatus == "case"

set.seed(1)
# ADAPT
library(ADAPT)

adapt_output <- adapt(otu_table=count_table, metadata=metadata, covar="iscase",
                      adjust=NULL, prevalence_cutoff=0.05, taxa_are_rows=T, boot=F)
# adapt_output_boot <- adapt(otu_table=count_table_subset, metadata=metadata, covar="iscase",
#                            adjust=NULL, prevalence_cutoff=0.05, genes_are_rows=T, boot=T, n_boot_gene=100,
#                            boot_replicate=1000, alpha=0.05)
adapt_DAtaxa <- adapt_output$DA_Taxa
adapt_details <- adapt_output$P_Value
adapt_summary <- data.frame(Taxa = adapt_DAtaxa, ADAPT=adapt_details[adapt_DAtaxa,  "effect"] > 0)
adapt_summary$ADAPT <- 2*(adapt_summary$ADAPT-0.5)
# adapt_reference_df <- adapt_result_df %>% filter(Gene %in% adapt_output_boot$Reference_Gene)
# adapt_DA_df <- adapt_result_df %>% filter(Gene %in% adapt_output_boot$DA_Gene)
adapt_summary$Taxa <- str_replace_all(adapt_summary$Taxa,"\\."," ")



# Aldex
library(ALDEx2)
aldex_result <- aldex(reads=count_table,
                         conditions=metadata$iscase)
aldex_DAtaxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]
aldex_summary <- data.frame(Taxa = aldex_DAtaxa, ALDEX=aldex_result[aldex_DAtaxa,  "effect"] > 0)
aldex_summary$ALDEX <- 2*(aldex_summary$ALDEX-0.5)
aldex_summary$Taxa <- str_replace_all(aldex_summary$Taxa,"\\."," ")


# Maaslin2
library(Maaslin2)
count_df <- data.frame(t(count_table))
maaslin2_output <-invisible(Maaslin2(input_data=count_df, input_metadata=metadata,
                                     output = file.path(folder, 'maaslin2_output'),
                                     min_prevalence=0.05,
                                     max_significance = 0.05,
                                     fixed_effects=c("iscase"),
                                     cores=4, plot_heatmap = F,plot_scatter = F))
maaslin2_result <- maaslin2_output$results
maaslin2_DAtaxa <- maaslin2_result$feature[maaslin2_result$qval < 0.05]
library(dplyr)
maaslin2_summary <- maaslin2_result %>% filter(qval < 0.05) %>%
  dplyr::select(feature, coef) %>% mutate(Taxa=feature, maaslin2=coef > 0, .keep="none")
maaslin2_summary$Taxa <-  str_replace_all(maaslin2_summary$Taxa,"\\."," ")
maaslin2_summary$maaslin2 <- 2*(maaslin2_summary$maaslin2 - 0.5)


# MetagenomeSeq, no DA taxa detected
library(metagenomeSeq)
MRmetadata <- AnnotatedDataFrame(metadata)
MRobj <- newMRexperiment(counts=count_table,
                         phenoData=MRmetadata)

MRobj <- wrenchNorm(MRobj, condition=metadata$CaseStatus)
mod <- model.matrix(~iscase, data=pData(MRobj))
MRfit <- fitFeatureModel(MRobj, mod)
MR_coefficients <- MRcoefs(MRfit, number=nrow(count_table), adjustMethod="BH")
MR_DAtaxa <- rownames(MR_coefficients)[MR_coefficients$adjPvalues < 0.05] |> na.omit()


# DACOMP, no DA taxa detected
library(dacomp)
set.seed(1)
selected_references <- dacomp.select_references(
  X = count_table, verbose = F, run_in_parallel = T, Nr.Cores=4,
  minimal_TA=10)
updated_selected_references <- dacomp.validate_references(X = t(count_table), #Counts
                                                          Y = metadata$iscase, #Traits
                                                          ref_obj = selected_references, #reference checked, must be object from dacomp.select_references(...)
                                                          test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                          Q_validation = 0.1, #FDR level for checking
                                                          Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                          Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                          NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                          Verbose = F)
dacomp_output <- dacomp.test(X = t(count_table), #counts data
                             y = metadata$iscase, #phenotype in y argument
                             # obtained from dacomp.select_references(...):
                             ind_reference_taxa = updated_selected_references,
                             test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                             disable_DSFDR = T, nr_perm=4000,
                             q=0.05,
                             verbose = F)

dacomp_DAtaxa <- taxa_names[dacomp_output$p.values.test.adjusted < 0.05] |>
  na.omit()



# ZicoSeq
library(GUniFrac)

zicoseq_output <- ZicoSeq(meta.dat=metadata,
                          feature.dat=count_table,
                          grp.name="iscase",
                          is.fwer=T,
                          perm.no=1999)
zicoseq_result <- cbind.data.frame(p.raw = zicoseq_output$p.raw,
                                   p.adj.fdr = zicoseq_output$p.adj.fdr,
                                   p.adj.fwer = zicoseq_output$p.adj.fwer,
                                   effect_sizes =zicoseq_output$coef.list[[1]][2,])
zicoseq_result$Taxa <- rownames(zicoseq_result)
zicoseq_summary <- zicoseq_result %>% filter(p.adj.fdr < 0.05) %>%
  mutate(Taxa=Taxa, zicoseq=effect_sizes > 0, .keep="none")
zicoseq_summary$zicoseq <- 2*(zicoseq_summary$zicoseq - 0.5)
zicoseq_summary$Taxa <- str_replace_all(zicoseq_summary$Taxa,"\\."," ")



#ANCOM
suppressMessages(library(ANCOMBC))
suppressMessages(library(phyloseq))
phyobj <- phyloseq(otu_table(count_table, taxa_are_rows = T),
                   sample_data(metadata))

ancom_output <- ancom(data=phyobj, p_adj_method="BH", prv_cut=0.05,
                      main_var="iscase", n_cl=4)
ancom_result <- ancom_output$res
ancom_avgcoef <- colMeans(ancom_output$beta_data)
ancom_summary <- data.frame(Taxa = ancom_output$res$taxon[ancom_output$res$detected_0.7])
ancom_summary$coef <- ancom_avgcoef[ancom_summary$Taxa]
ancom_summary <- ancom_summary %>% mutate(Taxa=Taxa, ANCOM=coef > 0, .keep="none")
ancom_summary$ANCOM <- 2*(ancom_summary$ANCOM - 0.5)
ancom_summary$Taxa <- str_replace_all(ancom_summary$Taxa,"\\."," ")





# ANCOMBC
ancombc_output <- ancombc2(phy_metag_plaque, fix_formula = 'iscase', p_adj_method='BH', global=T,
                           group='iscase', struc_zero=FALSE, prv_cut=0.05, n_cl=4)

ancombc_result <- ancombc_output$res
ancombc_summary <- ancombc_result %>% filter(diff_iscaseTRUE) %>%
  mutate(Taxa=taxon, ANCOMBC=lfc_iscaseTRUE > 0, .keep="none")
ancombc_summary$ANCOMBC <- 2*(ancombc_summary$ANCOMBC - 0.5)
ancombc_summary$Taxa <- str_replace_all(ancombc_summary$Taxa,"\\."," ")

# LinDA
library(LinDA)
linda_output <- linda(otu.tab=count_table,
                      meta=data.frame(metadata),
                      formula="~iscase",
                      n.cores=4)
linda_result <- linda_output$output$iscaseTRUE
linda_result$Taxa <- rownames(linda_result)
linda_summary <- linda_result %>% filter(padj < 0.05) %>%
  mutate(Taxa=Taxa, linda=log2FoldChange > 0, .keep="none")
linda_summary$linda <- 2*(linda_summary$linda - 0.5)
linda_summary$Taxa <- str_replace_all(linda_summary$Taxa,"\\."," ")

# aggregate outputs
library(purrr)
combined_result <- purrr::reduce(list(adapt_summary, aldex_summary, maaslin2_summary, zicoseq_summary, ancom_summary, ancombc_summary, linda_summary), 
                                 dplyr::full_join, by = 'Taxa')
combined_result$metagenomeSeq <- NA
combined_result$DACOMP <- NA
combined_result[is.na(combined_result)] <- 0

write.csv(combined_result, file.path(folder, "method_comparison_WGS_plaque.csv"),
          row.names=F)
write.csv(adapt_details, file.path(folder, "ADAPT_WGS_plaque_result.csv"),
          row.names=F)


