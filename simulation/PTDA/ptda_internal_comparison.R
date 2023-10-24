

rm(list=ls())
library(dplyr)
simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
polda_folder <- file.path(simulation_folder, 'POLDA')


mode <- "mix"

example_DA <- readRDS(file.path(polda_folder, "mix", "lognormal", 
                             "polda_n100_p10_d100000_t0.5_s20.rds"))
pval_DA_df <- example_DA$polda_result$P_Value
hist(pval_DA_df$pval, nclass=40)

example_null <- readRDS(file.path(polda_folder, "mix", "lognormal", 
                                  "polda_n100_p0_d100000_t0.5_s30.rds"))
pval_null_df <- example_null$polda_result$P_Value

hist(pval_null_df$pval, nclass=30)

combined_pvals <- data.frame(Pvals=c(pval_DA_df$pval, pval_null_df$pval),
                             Type=c(rep("With DA Taxa", 991), rep("Null", 988)))

ggplot(combined_pvals, aes(x=Pvals)) + geom_histogram(binwidth = 0.02, color="black", fill="white") + 
  facet_grid(cols=vars(Type)) + xlab("Raw Pvals") + ylab("Count") +
  theme_bw()+theme(text = element_text(size=12))


ID <- seq(1, 100)
tau <- c(0, 0.5, 1)
assumption <- c("lognormal", "loglogistic")

performance_df <- expand.grid(ID, tau, assumption)
colnames(performance_df) <- c("ID", "Penalty", "Assumption")
performance_df$reftaxa_error <- 0
performance_df$AUROC <- 0
performance_df$FDR <- 0
performance_df$Power <- 0

for (rowID in 1:nrow(performance_df)){
  tau <- performance_df$Penalty[rowID]
  seed <- performance_df$ID[rowID]
  assumption <- performance_df$Assumption[rowID]
  input_filename <- sprintf("polda_n100_p10_d20000_t%.1f_s%d.rds", tau, seed)
  model_summary <- readRDS(file.path(polda_folder, mode, assumption, input_filename))
  model_performance <- model_summary$performance
  performance_df$Power[rowID] <- model_performance$Power
  performance_df$reftaxa_error[rowID] <- model_performance$reftaxa_error
  performance_df$AUROC[rowID] <- model_performance$AUROC
  performance_df$FDR[rowID] <- model_performance$FDR
}

performance_df$Penalty <- sprintf("Tau=%.1f", performance_df$Penalty)
performance_df$Penalty <- factor(performance_df$Penalty, 
                                 levels=c("Tau=0.0", "Tau=0.5", "Tau=1.0"))



library(ggplot2)

FDR_plot <- ggplot(performance_df, aes(x=Assumption, y=FDR))  + geom_violin() +  geom_boxplot(width=0.2)+
  facet_grid(cols=vars(Penalty)) + geom_hline(yintercept=0.05, linetype="dashed", color="red", linewidth=1)+ylim(0, 0.2)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=14))

power_plot <- ggplot(performance_df, aes(x=Assumption, y=Power))  + geom_violin() +  geom_boxplot(width=0.2)+
  facet_grid(cols=vars(Penalty)) +ylim(0, 0.7)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=14))


library(cowplot)
plot_grid(FDR_plot, power_plot, align="v", ncol=1)



