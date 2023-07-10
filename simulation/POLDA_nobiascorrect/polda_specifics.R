# illustrate the details of POLDA using two examples
library(ggplot2)
library(cowplot)

rm(list=ls())

data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
polda_folder <- '/home/wangmk/MDAWG/POLDA/simulation/POLDA'

source("/home/wangmk/MDAWG/POLDA/POLDA/POLDA.R")
source(file.path(polda_folder, "utils.R"))


# The first example is a dataset with no DA taxa
ID <- 70
AGP_null <- readRDS(file.path(data_folder, "null", 
                              sprintf("AGP_simulation_null_%d.rds", ID)))
polda_null <- polda(otu_table=AGP_null$count_mat,
                    metadata=AGP_null$sample_metadata,
                    covar="X")
pval_null_df <- polda_null$P_Value

pval_hist_null <- ggplot(pval_null_df, aes(x=pval)) + 
  geom_histogram(binwidth=0.05, color="black", fill="white")+theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("P Values (No DA Taxa)") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))+
  xlim(-0.05, 1.05)

effect_estimate_null <- ggplot(pval_null_df, aes(x=effect))+
  geom_histogram(binwidth=0.2, color="black", fill="white")+
  geom_vline(xintercept=0, linetype="dashed", color="red")+ theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("Log Fold Change (No DA Taxa)") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))+
  xlim(-4, 4)


# The second example is a dataset with DA taxa
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                              sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))

countmat_unbalanced_rare_high <- AGP_unbalanced_rare_high$count_mat
sample_metadata_rare_high <- AGP_unbalanced_rare_high$sample_metadata
taxa_names <- rownames(countmat_unbalanced_rare_high)
reftaxa <- taxa_names

relabd_result_unbalanced_rare_high <- reference_GLM(count_data=countmat_unbalanced_rare_high,
                                                    metadata = sample_metadata_rare_high,
                                                    covar="X", reftaxa = reftaxa)

pval_hist_DA <- ggplot(relabd_result_unbalanced_rare_high, aes(x=pval))+ 
  geom_histogram(binwidth=0.05, color="black", fill="white")+theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("P Values (20% DA Taxa)") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))+
  xlim(-0.05, 1.05)

effect_estimate_DA <- ggplot(relabd_result_unbalanced_rare_high, aes(x=effect))+
  geom_histogram(binwidth=0.2, color="black", fill="white")+
  geom_vline(xintercept=0, linetype="dashed", color="red")+ theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("Log Fold Change (20% DA Taxa)") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))+
  xlim(-4, 4)

plot_grid(pval_hist_null, effect_estimate_null, pval_hist_DA,
          effect_estimate_DA, ncol=2)

## select a subset of taxa as potential reference set
estimated_effect <- relabd_result_unbalanced_rare_high$effect
names(estimated_effect) <- reftaxa
distance2med <- abs(estimated_effect - median(estimated_effect))
sorted_distance <- sort(distance2med)
ordered_taxanames <- names(sorted_distance)
reftaxa <- ordered_taxanames[1:(length(ordered_taxanames)/2)]
lower_bound <- min(estimated_effect[reftaxa])
upper_bound <- max(estimated_effect[reftaxa])

effect_estimate_DA + 
  geom_vline(xintercept=lower_bound, linetype="dashed", color="blue")+
  geom_vline(xintercept=upper_bound, linetype="dashed", color="blue")



relabd_result_subset <- reference_GLM(count_data=countmat_unbalanced_rare_high,
                                      metadata = sample_metadata_rare_high,
                                      covar="X", reftaxa = reftaxa)


pval_hist_DA_subset <- ggplot(relabd_result_subset, aes(x=pval))+ 
  geom_histogram(binwidth=0.05, color="black", fill="white")+theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("P Values Among Reference Taxa") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))+
  xlim(-0.05, 1.05)

effect_estimate_DA_subset <- ggplot(relabd_result_subset, aes(x=effect))+
  geom_histogram(binwidth=0.05, color="black", fill="white")+
  geom_vline(xintercept=0, linetype="dashed", color="red")+ theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("Log Fold Change Among Reference Taxa") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))

plot_grid(pval_hist_DA_subset, effect_estimate_DA_subset, ncol=2)


relabd_result_subset_complement <- reference_GLM(count_data=countmat_unbalanced_rare_high,
                                      metadata = sample_metadata_rare_high,
                                      covar="X", reftaxa = reftaxa, complement = TRUE)

relabd_result_final <- rbind(relabd_result_subset,
                             relabd_result_subset_complement)

pval_hist_DA_final <- ggplot(relabd_result_final, aes(x=pval))+ 
  geom_histogram(binwidth=0.05, color="black", fill="white")+theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("Final P Values") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))+
  xlim(-0.05, 1.05)

effect_estimate_DA_final <- ggplot(relabd_result_final, aes(x=effect))+
  geom_histogram(binwidth=0.2, color="black", fill="white")+
  geom_vline(xintercept=0, linetype="dashed", color="red")+ theme_bw()+
  xlab("") + ylab("Frequency") + ggtitle("Final Log Fold Change") +
  theme(text=element_text(size=11), plot.title=element_text(hjust=0.5))+
  xlim(-4, 4)

plot_grid(pval_hist_DA_final, effect_estimate_DA_final, ncol=2)


performance <- readRDS(file.path(polda_folder, "DA_rare", "high",
                                 sprintf("summary_unbalanced_rare_high_%d.rds", ID)))
