
rm(list=ls())
library(dispmod)
library(phyloseq)
library(survival)
folder <- '/home/wangmk/UM/Research/MDAWG/POLDA_example/real_data/ECC'
output_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA_example/real_data/ECC/16S'

load(file.path(folder, 'phyasv_meta.Rdata'))

meta12month <- data.frame(sample_data(ps_meta_nodups))
count_mat <- otu_table(ps_meta_nodups)@.Data
taxa_species <- tax_table(ps_meta_nodups)
cnames <- colnames(meta12month)

library(dplyr)

meta12month_subset <- meta12month %>% dplyr::select(FoxCavID, BabySubjectID, Visit, AgeAtExamMonths, Case, CaseStatus,
                                       Decay, IncidentVisit) %>%
  arrange(BabySubjectID, Visit)

meta12month <- meta12month_subset %>% mutate(visitto12 = abs(AgeAtExamMonths - 12)) %>%
  group_by(BabySubjectID) %>% filter(visitto12 == min(visitto12)) %>% filter(visitto12 <=2) %>%
  filter(Visit != IncidentVisit)
meta12month <- meta12month %>% group_by(BabySubjectID) %>% filter(row_number() == 1)
count_mat_12month <- count_mat[meta12month$FoxCavID, ]

prevalence <- colMeans(count_mat_12month != 0)

count_mat_12month <- count_mat_12month[, prevalence > 0.05]
# save a copy with zeros imputed with ones
count_mat_12month_copy <- count_mat_12month
count_mat_12month_copy[count_mat_12month_copy == 0] <- 1

meta12month$iscase <- meta12month$Case == "Case"
meta12month <- as.data.frame(meta12month)
rownames(meta12month) <- meta12month$FoxCavID

# POLDA
library(POLDA)
taxa_names <- colnames(count_mat_12month)
res.polda <- polda(otu_table =t(count_mat_12month), metadata=meta12month,
                      covar="iscase", prevalence_cutoff = 0.05, alpha=0.05,
                   ratio_model="loglogistic", firth=T)


refgene <- res.polda$Reference_Taxa
refcounts <- rowSums(count_mat_12month[, refgene])
polda_df <- res.polda$P_Value
polda_df$glm_effect <- 0
polda_df$glm_pval <- 0
polda_df$glm_padj <- 0
polda_df$cox_effect <- 0
polda_df$cox_pval <- 0
polda_df$cox_padj <- 0
for (j in 1:ncol(count_mat_12month)){
  
  singlegene_count <- count_mat_12month[, polda_df$Taxon[j]]
  singlegene_ratios <- rep(0, length(singlegene_count))
  exist <- singlegene_count > 0
  singlegene_ratios[exist] <- refcounts[exist] / singlegene_count[exist] 
  singlegene_ratios[!exist] <- refcounts[!exist]+1
  
  
  logit_model <- glm(cbind(singlegene_count, refcounts) ~ meta12month$iscase, family="binomial")
  overdispersed_logit_model <- glm.binomial.disp(logit_model, verbose=F) |> summary()
  polda_df$glm_effect[j] <- overdispersed_logit_model$coefficients[2, 1]
  polda_df$glm_pval[j] <- overdispersed_logit_model$coefficients[2, 4]
  
  cox_model <- coxph(Surv(singlegene_ratios, exist) ~ meta12month$iscase) |> summary()
  polda_df$cox_effect[j] <- cox_model$coefficients[1, 1]
  polda_df$cox_pval[j] <- cox_model$coefficients[1, 5]
  
}


polda_df$glm_padj <- p.adjust(polda_df$glm_pval, method="BH")
polda_df$cox_padj <- p.adjust(polda_df$cox_pval, method="BH")
polda_df$IsRef <- polda_df$Taxon %in% res.polda$Reference_Taxa
polda_df$IsDA <- polda_df$Taxon %in% res.polda$DA_taxa

res.polda$P_Value <- polda_df

write.csv(polda_df, file.path(output_folder, "polda_result.csv"), row.names=F)

loglogistic_DAtaxa <- res.polda$DA_taxa
cox_DAtaxa <- res.polda$P_Value$Taxon[which(res.polda$P_Value$cox_padj < 0.05)]
GLM_DAtaxa <- res.polda$P_Value$Taxon[which(res.polda$P_Value$glm_padj < 0.05)]

# LOCOM
library(LOCOM)
res.locom <- locom(otu.table = count_mat_12month, Y = meta12month$iscase, 
                   fdr.nominal = 0.05, prev.cut=0.05,
                   n.perm.max = 20000, n.rej.stop = 100, n.cores = 2)	

saveRDS(res.locom, file.path(output_folder, "LOCOM_16S.rds"))
locom_DAtaxa <- res.locom$detected.otu

# DACOMP
library(dacomp)
dacomp.ref.otus = dacomp.select_references(X = count_mat_12month, median_SD_threshold = 0, verbose = F)
# test.method = ifelse(length(unique(Y)) == 2, DACOMP.TEST.NAME.WILCOXON, DACOMP.TEST.NAME.SPEARMAN)
set.seed(2023)
res.dacomp <- dacomp.test(X = count_mat_12month, y = meta12month$iscase, nr_perm=5000, 
                           ind_reference_taxa = dacomp.ref.otus, 
                           test = DACOMP.TEST.NAME.WILCOXON, verbose = F, q = 0.05)
dacomp_DAtaxa <- taxa_names[res.dacomp$dsfdr_rejected]

saveRDS(res.dacomp, file.path(output_folder, "DACOMP_16S.rds"))



# ALDEx2
library(ALDEx2)
res.aldex <- aldex(reads=t(count_mat_12month),
                  meta12month$iscase, denom="iqlr")
aldex_DAtaxa <- taxa_names[res.aldex$wi.eBH < 0.05]
saveRDS(res.aldex, file.path(output_folder, "aldex_16S.rds"))


# ANCOM-BC
library(ANCOMBC)
otu.phylo <- otu_table(count_mat_12month, taxa_are_rows = FALSE)
meta.phylo <- sample_data(meta12month)
# tax.phylo <- tax_table(tax.table)
phylo.obj <- phyloseq(otu.phylo, meta.phylo)

res.ancombc <- ancombc(phyloseq = phylo.obj, formula = "iscase",
                       p_adj_method = "BH", zero_cut = 0.95, lib_cut = 1000, # 0.2: LOCOM filter; 0.1: ANCOM filter
                       group = "iscase", struc_zero = FALSE, neg_lb = FALSE,
                       tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)
ancombc_DAtaxa <- taxa_names[res.ancombc$res$q_val < 0.05]

ANCOMBC_taxa <- which(res.ancombc$res$q_val$iscaseTRUE < 0.05)
saveRDS(res.ancombc, file.path(output_folder, "ANCOMBC_16S.rds"))

DAresults <- list(polda_loglogistic=loglogistic_DAtaxa,
                  polda_GLM=GLM_DAtaxa,
                  polda_cox=cox_DAtaxa,
                  ALDeX2=aldex_DAtaxa,
                  ANCOMBC=ancombc_DAtaxa,
                  DACOMP=dacomp_DAtaxa,
                  LOCOM=locom_DAtaxa)

saveRDS(DAresults, file.path(output_folder, "DA_results.rds"))

count_intersect <- matrix(0, nrow=7, ncol=7)
colnames(count_intersect) <- c("loglogistic", "GLM", "Cox", "ALDEx2", "ANCOMBC", "DACOMP", "LOCOM")
rownames(count_intersect) <- c("loglogistic", "GLM", "Cox", "ALDEx2", "ANCOMBC", "DACOMP", "LOCOM")
for (j in 1:6){
  for(k in (j+1):7){
    count_intersect[j, k] <- length(intersect(DAresults[[j]], DAresults[[k]]))
  }
}
for (j in 1:7){
  count_intersect[j, j] <- length(DAresults[[j]])
}
write.csv(count_intersect, file.path(output_folder, "DA_result_comparison.csv"))

