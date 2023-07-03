rm(list=ls())
library(locfdr)

POLDA_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/POLDA"
LOCOM_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/LOCOM"
metagenomeseq_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/MetagenomeSeq"
data_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/data"

polda_model <- readRDS(file.path(POLDA_folder, "DA_rare", "medium",
                                  "summary_unbalanced_rare_medium_1.rds"))

locom_model <- readRDS(file.path(LOCOM_folder, "DA_rare", "medium",
                                 "summary_unbalanced_rare_medium_1.rds"))


metagenomeseq_model <- readRDS(file.path(metagenomeseq_folder, "DA_rare", "medium",
                                         "summary_unbalanced_rare_medium_1.rds"))


taxa_names <- sprintf("Taxon_%d", seq(1, 1000))
polda_result <- polda_model$polda_result
pval_df <- polda_result$P_Value
polda_raw_pvals <- polda_result$P_Value$pval
names(polda_raw_pvals) <- pval_df$Taxon
polda_raw_pvals <- polda_raw_pvals[taxa_names]
locom_raw_pvals <- as.vector(locom_model$locom_result$p.otu)
rawpval_df <- data.frame(polda=polda_raw_pvals,
                         locom=locom_raw_pvals)
hist(locom_raw_pvals, nclass=20)

adjusted_locom_pvals <- p.adjust(locom_raw_pvals, method="BH")

filter <- !(polda_result$P_Value$Taxon %in% polda_result$Reference_Taxa)

sum(p.adjust(raw_pvals[filter], method="BH") < 0.05)

simulated_data <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                    "AGP_simulation_unbalanced_rare_medium_1.rds"))


