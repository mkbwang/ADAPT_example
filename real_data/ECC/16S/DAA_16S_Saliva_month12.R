library(dplyr)
library(stringr)
folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/16S"

phy_asv_month12 <- readRDS(file.path(folder, "phyloseq", "phyasv_visit12.rds"))

phy_asv_month12 <- phy_asv_month12 %>% filter_taxa(function(ct) mean(ct > 0) > 0.05, TRUE)
count_table <- otu_table(phy_asv_month12)@.Data
taxa_names <- colnames(count_table)
metadata <- data.frame(sample_data(phy_asv_month12))
metadata$iscase <- 1*(metadata$CaseEver == "Case")
metadata$female <- 1*(metadata$host_sex == "female")
metadata$WV <- 1*(metadata$geo_loc_name == "USA: WV")


set.seed(1)
# ADAPT
library(ADAPT)
adapt_output <- adapt(input_data=phy_asv_month12, cond.var = "CaseEver", base.cond = "Control")
adapt_DAtaxa <- adapt_output@signal
adapt_details <- adapt_output@details
adapt_summary <- data.frame(Taxa = adapt_DAtaxa, ADAPT=adapt_details[adapt_DAtaxa,  "log10foldchange"] > 0)
adapt_summary$ADAPT <- 2*(adapt_summary$ADAPT-0.5)



# Aldex
library(ALDEx2)
aldex_result <- aldex(reads=t(count_table),
                      conditions=metadata$iscase)
aldex_DAtaxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]
aldex_summary <- data.frame(Taxa = aldex_DAtaxa, ALDEx2=aldex_result[aldex_DAtaxa,  "effect"] > 0)
aldex_summary$ALDEx2 <- 2*(aldex_summary$ALDEx2-0.5)


# Maaslin2
library(Maaslin2)
count_df <- data.frame(count_table)
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
  dplyr::select(feature, coef) %>% mutate(Taxa=feature, Maaslin2=coef > 0, .keep="none")
maaslin2_summary$Maaslin2 <- 2*(maaslin2_summary$Maaslin2 - 0.5)


# MetagenomeSeq, no DA taxa detected
library(metagenomeSeq)
MRmetadata <- AnnotatedDataFrame(metadata)
MRobj <- newMRexperiment(counts=t(count_table),
                         phenoData=MRmetadata)

MRobj <- wrenchNorm(MRobj, condition=metadata$CaseEver)
mod <- model.matrix(~iscase, data=pData(MRobj))
MRfit <- fitFeatureModel(MRobj, mod)
MR_coefficients <- MRcoefs(MRfit, number=ncol(count_table), adjustMethod="BH")
MR_coefficients$Taxa <- rownames(MR_coefficients)
MR_DAtaxa <- rownames(MR_coefficients)[MR_coefficients$adjPvalues < 0.05] |> na.omit()
MR_summary <- MR_coefficients %>% filter(adjPvalues < 0.05) %>%
  mutate(Taxa=Taxa, metagenomeSeq=logFC > 0, .keep="none")
MR_summary$metagenomeSeq <- 2*(MR_summary$metagenomeSeq - 0.5)


# DACOMP, no DA taxa detected
library(dacomp)
set.seed(1)
selected_references <- dacomp.select_references(
  X = count_table, verbose = F, run_in_parallel = T, Nr.Cores=4,
  minimal_TA=10)
updated_selected_references <- dacomp.validate_references(X = count_table, #Counts
                                                          Y = metadata$iscase, #Traits
                                                          ref_obj = selected_references, #reference checked, must be object from dacomp.select_references(...)
                                                          test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                          Q_validation = 0.1, #FDR level for checking
                                                          Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                          Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                          NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                          Verbose = F)
dacomp_output <- dacomp.test(X = count_table, #counts data
                             y = metadata$iscase, #phenotype in y argument
                             # obtained from dacomp.select_references(...):
                             ind_reference_taxa = updated_selected_references,
                             test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                             disable_DSFDR = T, nr_perm=4000,
                             q=0.05,
                             verbose = F)
DA_filter <- dacomp_output$p.values.test.adjusted < 0.05 & !is.na(dacomp_output$p.values.test.adjusted)
Direction <- dacomp_output$effect_size_estimates[DA_filter]
dacomp_summary <- data.frame(Taxa = taxa_names[DA_filter])
dacomp_summary$DACOMP = grepl("0<=1", Direction)
dacomp_summary$DACOMP = 2*(dacomp_summary$DACOMP - 0.5)


# ZicoSeq
library(GUniFrac)

zicoseq_output <- ZicoSeq(meta.dat=metadata,
                          feature.dat=t(count_table),
                          grp.name="iscase",
                          is.fwer=T,
                          perm.no=1999)
zicoseq_result <- cbind.data.frame(p.raw = zicoseq_output$p.raw,
                                   p.adj.fdr = zicoseq_output$p.adj.fdr,
                                   p.adj.fwer = zicoseq_output$p.adj.fwer,
                                   effect_sizes =zicoseq_output$coef.list[[1]][2,])
zicoseq_result$Taxa <- rownames(zicoseq_result)
zicoseq_summary <- zicoseq_result %>% filter(p.adj.fdr < 0.05) %>%
  mutate(Taxa=Taxa, ZicoSeq=effect_sizes > 0, .keep="none")
zicoseq_summary$ZicoSeq <- 2*(zicoseq_summary$ZicoSeq - 0.5)



#ANCOM
suppressMessages(library(ANCOMBC))
suppressMessages(library(phyloseq))
begin <- proc.time()
ancom_output <- ancom(data=phy_asv_month12, p_adj_method="BH", prv_cut=0.05,
                      main_var="CaseEver", n_cl=4)
duration <- proc.time() - begin
ancom_result <- ancom_output$res
ancom_avgcoef <- rowMeans(ancom_output$beta_data)
ancom_summary <- data.frame(Taxa = ancom_output$res$taxon[ancom_output$res$detected_0.7])
ancom_summary$coef <- ancom_avgcoef[ancom_summary$Taxa]
ancom_summary <- ancom_summary %>% mutate(Taxa=Taxa, ANCOM=coef > 0, .keep="none")
ancom_summary$ANCOM <- 2*(ancom_summary$ANCOM - 0.5)


# ANCOMBC
ancombc_output <- ancombc2(phy_asv_month12, fix_formula = 'CaseEver', p_adj_method='BH', global=T,
                           group='CaseEver', struc_zero=FALSE, prv_cut=0.05, n_cl=4)

ancombc_result <- ancombc_output$res
ancombc_summary <- ancombc_result %>% filter(diff_CaseEverControl) %>%
  mutate(Taxa=taxon, ANCOMBC=lfc_CaseEverControl < 0, .keep="none")
ancombc_summary$ANCOMBC <- 2*(ancombc_summary$ANCOMBC - 0.5)


# LinDA
suppressMessages(library(MicrobiomeStat))
linda_output <- linda(feature.dat=t(count_table),
                      meta.dat=metadata,
                      formula="~iscase",
                      n.cores=4)
linda_result <- linda_output$output$iscase
linda_result$Taxa <- rownames(linda_result)
linda_summary <- linda_result %>% filter(padj < 0.05) %>%
  mutate(Taxa=Taxa, LinDA=log2FoldChange > 0, .keep="none")
linda_summary$LinDA <- 2*(linda_summary$LinDA - 0.5)


# aggregate outputs
library(purrr)
combined_result <- purrr::reduce(list(adapt_summary, aldex_summary, maaslin2_summary, MR_summary, dacomp_summary,
                                      zicoseq_summary, ancom_summary, ancombc_summary, linda_summary), 
                                 dplyr::full_join, by = 'Taxa')
combined_result[is.na(combined_result)] <- 0

taxonomy_table <- tax_table(phy_asv_month12)@.Data
selected_taxonomies <- taxonomy_table[combined_result$Taxa, ] |> as.data.frame()

all_DNAstrings <- as.character(phy_asv_month12@refseq)
selected_DNAstrings <- all_DNAstrings[combined_result$Taxa]
selected_taxonomies$DNAstrings <- selected_DNAstrings
combined_result <- cbind(combined_result, selected_taxonomies)


write.csv(combined_result, file.path(folder, "method_comparison_16S_saliva_12month.csv"),
          row.names=F)
write.csv(adapt_details, file.path(folder, "ADAPT_16S_saliva_12month_result.csv"),
          row.names=F)


