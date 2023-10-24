

rm(list=ls())
library(dplyr)
simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
dacomp_folder <- file.path(simulation_folder, 'DACOMP')



dacomp_performance <- data.frame(ID=seq(1, 100),
                                reftaxa_error=rep(0, 100),
                                AUROC=rep(0, 100),
                                FDR=rep(0, 100),
                                Power=rep(0, 100))


input_folder <- "mix"
# ID <- 30
# input_filename <- sprintf("dacomp_n100_p10_d100000_s%d.rds", ID)
# data_filename <- sprintf("sim_n100_p10_d100000_s%d.rds", ID)
# 
# example_output <- readRDS(file.path(dacomp_folder, input_folder, input_filename))
# dacomp_output <- example_output$dacomp_result
# 
# pval_df <- dacomp_output$P_Value
# pval_df_significant <- pval_df %>% filter(adjusted_pval < 0.05)
# 
# original_data <- readRDS(file.path(data_folder, input_folder, data_filename))
# taxa_info <- origi
# zero_proportions1 <- rowMeans(count_mat[FP_taxa, 1:50] == 0)
# zero_proportions2 <- rowMeans(count_mat[FP_taxa, 51:100] == 0)nal_data$taxa_info
# count_mat <- original_data$count_mat
# 
# true_DAtaxa <- rownames(original_data$taxa_info)[original_data$taxa_info$Logfold != 0]
# 
# fp_df <- pval_df_significant %>% filter(!Taxon %in% true_DAtaxa)
# 
# FP_taxa <- fp_df$Taxon
# FP_counts <- count_mat[FP_taxa, ]

for (ID in 1:100){
  input_filename <- sprintf("dacomp_n100_p10_d100000_s%d.rds", ID)
  model_summary <- readRDS(file.path(dacomp_folder, input_folder, input_filename))
  model_performance <- model_summary$performance
  dacomp_performance$Power[ID] <- model_performance$Power
  dacomp_performance$reftaxa_error[ID] <- model_performance$reftaxa_error
  dacomp_performance$AUROC[ID] <- model_performance$AUROC
  dacomp_performance$FDR[ID] <- model_performance$FDR
}

boxplot(dacomp_performance$FDR)
abline(h=0.05, col="blue")
