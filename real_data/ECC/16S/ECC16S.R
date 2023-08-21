
rm(list=ls())
library(dispmod)
library(survival)
folder <- '/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/processed'
output_folder <- '/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data/ECC/16S'

load(file.path(folder, 'phyasv_meta.Rdata'))

meta12month <- data.frame(sample_data(ps_meta_nodups))
count_mat <- otu_table(ps_meta_nodups)@.Data
cnames <- colnames(meta12month)

library(dplyr)

meta12month_subset <- meta12month %>% dplyr::select(FoxCavID, BabySubjectID, Visit, AgeAtExamMonths, Case, CaseStatus,
                                       Decay, IncidentVisit) %>%
  arrange(BabySubjectID, Visit)

meta12month <- meta12month_subset %>% mutate(visitto12 = abs(AgeAtExamMonths - 12)) %>%
  group_by(BabySubjectID) %>% filter(visitto12 == min(visitto12)) %>% filter(visitto12 <=2) %>%
  filter(Visit != IncidentVisit)
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

res.polda$P_Value <- polda_df

saveRDS(res.polda, file.path(output_folder, "polda_16S.rds"))




# LOCOM
library(LOCOM)
res.locom <- locom(otu.table = count_mat_12month, Y = meta12month$iscase, 
                   fdr.nominal = 0.05, prev.cut=0.05,
                   n.perm.max = 20000, n.rej.stop = 100, n.cores = 2)	

saveRDS(res.locom, file.path(output_folder, "LOCOM_16S.rds"))


# DACOMP
library(dacomp)
dacomp.ref.otus = dacomp.select_references(X = count_mat_12month, median_SD_threshold = 0, verbose = F)
# test.method = ifelse(length(unique(Y)) == 2, DACOMP.TEST.NAME.WILCOXON, DACOMP.TEST.NAME.SPEARMAN)
res.dacomp <- dacomp.test(X = count_mat_12month, y = meta12month$iscase, 
                           ind_reference_taxa = dacomp.ref.otus, 
                           test = DACOMP.TEST.NAME.WILCOXON, verbose = F, q = 0.05)

saveRDS(res.dacomp, file.path(output_folder, "DACOMP_16S.rds"))



# ALDEx2
library(ALDEx2)
res.aldex <- aldex(reads=t(count_mat_12month),
                  meta12month$iscase, denom="iqlr")
saveRDS(res.aldex, file.path(output_folder, "aldex_16S.rds"))


# ANCOM-BC
library(ANCOMBC)
otu.phylo <- otu_table(count_mat_12month, taxa_are_rows = FALSE)
meta.phylo <- sample_data(meta12month)
# tax.phylo <- tax_table(tax.table)
phylo.obj <- phyloseq(otu.phylo, meta.phylo)

res.ancombc <- ancombc(phyloseq = phylo.obj, formula = "iscase",
                       p_adj_method = "BH", prv_cut = 0, lib_cut = 1000, # 0.2: LOCOM filter; 0.1: ANCOM filter
                       group = "iscase", struc_zero = FALSE, neg_lb = FALSE,
                       tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

ANCOMBC_taxa <- which(res.ancombc$res$q_val$iscaseTRUE < 0.05)
saveRDS(res.ancombc, file.path(output_folder, "ANCOMBC_16S.rds"))

