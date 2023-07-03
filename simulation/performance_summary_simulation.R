rm(list=ls())
# performance_summary
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
data_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA/simulation/data'
# aldex2_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA/simulation/Aldex2'
ancombc_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA/simulation/ANCOMBC'
# DACOMP_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA/simulation/DACOMP'
POLDA_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/POLDA"
# LOCOM_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/LOCOM"
MetagenomeSeq_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA/simulation/MetagenomeSeq'
# deseq2_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/Deseq2"
# corncob_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/corncob"
# maaslin2_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/Maaslin2"
# RAIDA_folder <- "/home/wangmk/UM/Research/MDAWG/POLDA/simulation/RAIDA"

template_df <- data.frame(FDR=rep(0, 100),
                          Power=rep(0, 100),
                          AUC=rep(0, 100))

polda_performance <- list(low=template_df,
                          medium=template_df,
                          high=template_df)
# polda result summary
for (j in 1:100){

  DA_low <- readRDS(file.path(POLDA_folder, 
                                       "DA_mixed", 
                                       "low", 
                                       sprintf("summary_unbalanced_mixed_low_%d.rds", j)))
  polda_performance$low$FDR[j] <- DA_low$performance$FDR
  polda_performance$low$Power[j] <- DA_low$performance$Power
  polda_performance$low$AUC[j] <- DA_low$performance$AUROC
  
  DA_medium <- readRDS(file.path(POLDA_folder, 
                                 "DA_mixed", 
                                 "medium", 
                                 sprintf("summary_unbalanced_mixed_medium_%d.rds", j)))
  polda_performance$medium$FDR[j] <- DA_medium$performance$FDR
  polda_performance$medium$Power[j] <- DA_medium$performance$Power
  polda_performance$medium$AUC[j] <- DA_medium$performance$AUROC
  
  DA_high <- readRDS(file.path(POLDA_folder, 
                               "DA_mixed", 
                               "high", 
                               sprintf("summary_unbalanced_mixed_high_%d.rds", j)))
  polda_performance$high$FDR[j] <- DA_high$performance$FDR
  polda_performance$high$Power[j] <- DA_high$performance$Power
  polda_performance$high$AUC[j] <- DA_high$performance$AUROC
  
}

polda_performance$low$DA_prop <- "5%"
polda_performance$medium$DA_prop <- "10%"
polda_performance$high$DA_prop <- "20%"
polda_combined_performance <- rbind(polda_performance$low,
                                    polda_performance$medium,
                                    polda_performance$high)
polda_combined_performance$Method <- "POLDA"
polda_combined_performance$DA_prop <- factor(polda_combined_performance$DA_prop,
                                             levels=c("5%", "10%", "20%"))



# wrench + MSeq
MetagenomeSeq_performance <- list(low=template_df,
                                  medium=template_df,
                                  high=template_df)
# MetagenomeSeq result summary
for (j in 1:100){
  
  DA_low <- readRDS(file.path(MetagenomeSeq_folder, 
                              "DA_mixed", 
                              "low", 
                              sprintf("summary_unbalanced_mixed_low_%d.rds", j)))
  MetagenomeSeq_performance$low$FDR[j] <- DA_low$performance$FDR
  MetagenomeSeq_performance$low$Power[j] <- DA_low$performance$Power
  MetagenomeSeq_performance$low$AUC[j] <- DA_low$performance$AUROC
  
  DA_medium <- readRDS(file.path(MetagenomeSeq_folder, 
                                 "DA_mixed", 
                                 "medium", 
                                 sprintf("summary_unbalanced_mixed_medium_%d.rds", j)))
  MetagenomeSeq_performance$medium$FDR[j] <- DA_medium$performance$FDR
  MetagenomeSeq_performance$medium$Power[j] <- DA_medium$performance$Power
  MetagenomeSeq_performance$medium$AUC[j] <- DA_medium$performance$AUROC
  
  DA_high <- readRDS(file.path(MetagenomeSeq_folder, 
                               "DA_mixed", 
                               "high", 
                               sprintf("summary_unbalanced_mixed_high_%d.rds", j)))
  MetagenomeSeq_performance$high$FDR[j] <- DA_high$performance$FDR
  MetagenomeSeq_performance$high$Power[j] <- DA_high$performance$Power
  MetagenomeSeq_performance$high$AUC[j] <- DA_high$performance$AUROC
  
}

MetagenomeSeq_performance$low$DA_prop <- "5%"
MetagenomeSeq_performance$medium$DA_prop <- "10%"
MetagenomeSeq_performance$high$DA_prop <- "20%"
MetagenomeSeq_combined_performance <- rbind(MetagenomeSeq_performance$low,
                                            MetagenomeSeq_performance$medium,
                                            MetagenomeSeq_performance$high)
MetagenomeSeq_combined_performance$Method <- "MetagenomeSeq"
MetagenomeSeq_combined_performance$DA_prop <- factor(MetagenomeSeq_combined_performance$DA_prop,
                                             levels=c("5%", "10%", "20%"))



ancombc_performance <- list(low=template_df,
                            medium=template_df,
                            high=template_df)
# ancombc result summary
for (j in 1:100){
  
  DA_low <- readRDS(file.path(ancombc_folder, 
                              "DA_mixed", 
                              "low", 
                              sprintf("summary_unbalanced_mixed_low_%d.rds", j)))
  ancombc_performance$low$FDR[j] <- DA_low$performance$FDR
  ancombc_performance$low$Power[j] <- DA_low$performance$Power
  ancombc_performance$low$AUC[j] <- DA_low$performance$AUROC
  
  DA_medium <- readRDS(file.path(ancombc_folder, 
                                 "DA_mixed", 
                                 "medium", 
                                 sprintf("summary_unbalanced_mixed_medium_%d.rds", j)))
  ancombc_performance$medium$FDR[j] <- DA_medium$performance$FDR
  ancombc_performance$medium$Power[j] <- DA_medium$performance$Power
  ancombc_performance$medium$AUC[j] <- DA_medium$performance$AUROC
  
  DA_high <- readRDS(file.path(ancombc_folder, 
                               "DA_mixed", 
                               "high", 
                               sprintf("summary_unbalanced_mixed_high_%d.rds", j)))
  ancombc_performance$high$FDR[j] <- DA_high$performance$FDR
  ancombc_performance$high$Power[j] <- DA_high$performance$Power
  ancombc_performance$high$AUC[j] <- DA_high$performance$AUROC
  
}

ancombc_performance$low$DA_prop <- "5%"
ancombc_performance$medium$DA_prop <- "10%"
ancombc_performance$high$DA_prop <- "20%"
ancombc_combined_performance <- rbind(ancombc_performance$low,
                                      ancombc_performance$medium,
                                      ancombc_performance$high)
ancombc_combined_performance$Method <- "ANCOMBC"
ancombc_combined_performance$DA_prop <- factor(ancombc_combined_performance$DA_prop,
                                                     levels=c("5%", "10%", "20%"))

combined_performance <- rbind(polda_combined_performance, ancombc_combined_performance,
                              MetagenomeSeq_combined_performance)


library(ggplot2)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


FDR_plot <- ggplot(combined_performance, aes(x=DA_prop, y=FDR, color=Method)) + 
  geom_boxplot() + 
  scale_color_manual(values=cbp1)+
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  xlab("Proportion of DA Taxa") + 
  scale_y_continuous(breaks=seq(0, 0.3, 0.05), limits=c(0, 0.3))+
  ggtitle("False Discovery Rate")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5),
        axis.title.y=element_blank(),
        legend.position='bottom')


Power_plot <- ggplot(combined_performance, aes(x=DA_prop, y=Power, color=Method)) + 
  geom_boxplot() + 
  scale_color_manual(values=cbp1)+
  xlab("Proportion of DA Taxa") + 
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1))+
  ggtitle("Power")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5),
        axis.title.y=element_blank(),
        legend.position='bottom')




library(patchwork)
library(ggpubr)

ggarrange(FDR_plot, Power_plot, 
          ncol=2, common.legend=TRUE, 
          legend='bottom')


deseq2_performance <- list(FPR_null = rep(0, 100),
                           rare_5=template_df,
                           rare_10=template_df,
                           rare_20=template_df,
                           abundant_5=template_df,
                           abundant_10=template_df,
                           abundant_20=template_df)
# deseq2 result summary
for (j in 1:100){
  
  null_example <- readRDS(file.path(deseq2_folder, 
                                    "null", 
                                    sprintf("summary_null_%d.rds", j)))
  deseq2_performance$FPR_null[j] <- null_example$performance$FPR
  
  DA_rare_5 <- readRDS(file.path(deseq2_folder, 
                                 "DA_rare", 
                                 "low", 
                                 sprintf("summary_unbalanced_rare_low_%d.rds", j)))
  deseq2_performance$rare_5$FDR[j] <- DA_rare_5$performance$FDR
  deseq2_performance$rare_5$Power[j] <- DA_rare_5$performance$Power
  deseq2_performance$rare_5$AUC[j] <- DA_rare_5$performance$AUROC
  
  DA_rare_10 <- readRDS(file.path(deseq2_folder, 
                                  "DA_rare", 
                                  "medium", 
                                  sprintf("summary_unbalanced_rare_medium_%d.rds", j)))
  deseq2_performance$rare_10$FDR[j] <- DA_rare_10$performance$FDR
  deseq2_performance$rare_10$Power[j] <- DA_rare_10$performance$Power
  deseq2_performance$rare_10$AUC[j] <- DA_rare_10$performance$AUROC
  
  DA_rare_20 <- readRDS(file.path(deseq2_folder, 
                                  "DA_rare", 
                                  "high", 
                                  sprintf("summary_unbalanced_rare_high_%d.rds", j)))
  deseq2_performance$rare_20$FDR[j] <- DA_rare_20$performance$FDR
  deseq2_performance$rare_20$Power[j] <- DA_rare_20$performance$Power
  deseq2_performance$rare_20$AUC[j] <- DA_rare_20$performance$AUROC
  
  
  DA_abundant_5 <- readRDS(file.path(deseq2_folder, 
                                     "DA_abundant", 
                                     "low", 
                                     sprintf("summary_unbalanced_abundant_low_%d.rds", j)))
  deseq2_performance$abundant_5$FDR[j] <- DA_abundant_5$performance$FDR
  deseq2_performance$abundant_5$Power[j] <- DA_abundant_5$performance$Power
  deseq2_performance$abundant_5$AUC[j] <- DA_abundant_5$performance$AUROC
  
  DA_abundant_10 <- readRDS(file.path(deseq2_folder, 
                                      "DA_abundant", 
                                      "medium", 
                                      sprintf("summary_unbalanced_abundant_medium_%d.rds", j)))
  deseq2_performance$abundant_10$FDR[j] <- DA_abundant_10$performance$FDR
  deseq2_performance$abundant_10$Power[j] <- DA_abundant_10$performance$Power
  deseq2_performance$abundant_10$AUC[j] <- DA_abundant_10$performance$AUROC
  
  DA_abundant_20 <- readRDS(file.path(deseq2_folder, 
                                      "DA_abundant", 
                                      "high", 
                                      sprintf("summary_unbalanced_abundant_high_%d.rds", j)))
  deseq2_performance$abundant_20$FDR[j] <- DA_abundant_20$performance$FDR
  deseq2_performance$abundant_20$Power[j] <- DA_abundant_20$performance$Power
  deseq2_performance$abundant_20$AUC[j] <- DA_abundant_20$performance$AUROC
  
}
deseq2_performance$Method <- 'deseq2'



# Maaslin2
maaslin2_performance <- list(FPR_null = rep(0, 100),
                             rare_5=template_df,
                             rare_10=template_df,
                             rare_20=template_df,
                             abundant_5=template_df,
                             abundant_10=template_df,
                             abundant_20=template_df)
# maaslin2 result summary
for (j in 1:100){
  
  null_example <- readRDS(file.path(maaslin2_folder, 
                                    "null", 
                                    sprintf("summary_null_%d.rds", j)))
  maaslin2_performance$FPR_null[j] <- null_example$performance$FPR
  
  DA_rare_5 <- readRDS(file.path(maaslin2_folder, 
                                 "DA_rare", 
                                 "low", 
                                 sprintf("summary_unbalanced_rare_low_%d.rds", j)))
  maaslin2_performance$rare_5$FDR[j] <- DA_rare_5$performance$FDR
  maaslin2_performance$rare_5$Power[j] <- DA_rare_5$performance$Power
  maaslin2_performance$rare_5$AUC[j] <- DA_rare_5$performance$AUROC
  
  DA_rare_10 <- readRDS(file.path(maaslin2_folder, 
                                  "DA_rare", 
                                  "medium", 
                                  sprintf("summary_unbalanced_rare_medium_%d.rds", j)))
  maaslin2_performance$rare_10$FDR[j] <- DA_rare_10$performance$FDR
  maaslin2_performance$rare_10$Power[j] <- DA_rare_10$performance$Power
  maaslin2_performance$rare_10$AUC[j] <- DA_rare_10$performance$AUROC
  
  DA_rare_20 <- readRDS(file.path(maaslin2_folder, 
                                  "DA_rare", 
                                  "high", 
                                  sprintf("summary_unbalanced_rare_high_%d.rds", j)))
  maaslin2_performance$rare_20$FDR[j] <- DA_rare_20$performance$FDR
  maaslin2_performance$rare_20$Power[j] <- DA_rare_20$performance$Power
  maaslin2_performance$rare_20$AUC[j] <- DA_rare_20$performance$AUROC
  
  
  DA_abundant_5 <- readRDS(file.path(maaslin2_folder, 
                                     "DA_abundant", 
                                     "low", 
                                     sprintf("summary_unbalanced_abundant_low_%d.rds", j)))
  maaslin2_performance$abundant_5$FDR[j] <- DA_abundant_5$performance$FDR
  maaslin2_performance$abundant_5$Power[j] <- DA_abundant_5$performance$Power
  maaslin2_performance$abundant_5$AUC[j] <- DA_abundant_5$performance$AUROC
  
  DA_abundant_10 <- readRDS(file.path(maaslin2_folder, 
                                      "DA_abundant", 
                                      "medium", 
                                      sprintf("summary_unbalanced_abundant_medium_%d.rds", j)))
  maaslin2_performance$abundant_10$FDR[j] <- DA_abundant_10$performance$FDR
  maaslin2_performance$abundant_10$Power[j] <- DA_abundant_10$performance$Power
  maaslin2_performance$abundant_10$AUC[j] <- DA_abundant_10$performance$AUROC
  
  DA_abundant_20 <- readRDS(file.path(maaslin2_folder, 
                                      "DA_abundant", 
                                      "high", 
                                      sprintf("summary_unbalanced_abundant_high_%d.rds", j)))
  maaslin2_performance$abundant_20$FDR[j] <- DA_abundant_20$performance$FDR
  maaslin2_performance$abundant_20$Power[j] <- DA_abundant_20$performance$Power
  maaslin2_performance$abundant_20$AUC[j] <- DA_abundant_20$performance$AUROC
  
}
maaslin2_performance$Method <- 'maaslin2'



combined_results <- rbind(polda_performance,
                          MetagenomeSeq_performance,
                          deseq2_performance,
                          maaslin2_performance,
                          ancombc_performance)

FPR_results <- combined_results[, c("Method", "FPR_null")]
FDR_rare_results <- combined_results[, c("Method", "FDR_rare")]
Power_rare_results <- combined_results[, c("Method", "Power_rare")]
FDR_abundant_results <- combined_results[, c("Method", "FDR_abundant")]
Power_abundant_results <- combined_results[, c("Method", "Power_abundant")]

library(dplyr)
FPR_avg <- FPR_results %>% group_by(Method) %>% 
  summarise(FPR_avg = mean(FPR_null)) %>% arrange(desc(FPR_avg))
FPR_results$Method <- factor(FPR_results$Method,
                                levels=FPR_avg$Method)
FPR_comparison_plot <- ggplot(FPR_results, aes(x=Method, y=FPR_null)) + geom_violin() +
  geom_boxplot(width=0.02) + ylim(0, 0.2) +
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  xlab("") + ylab("False Positive Rate") +ggtitle("No DA Taxa")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5))+ coord_flip()


FDR_rare_avg <- FDR_rare_results %>% group_by(Method) %>% 
  summarise(FDR_rare_avg = mean(FDR_rare)) %>% arrange(desc(FDR_rare_avg))
FDR_rare_results$Method <- factor(FDR_rare_results$Method,
                             levels=FDR_rare_avg$Method)
FDR_rare_plot <- ggplot(FDR_rare_results, aes(x=Method, y=FDR_rare)) + geom_violin(scale="count") +
  geom_boxplot(width=0.05) + ylim(0, 0.5) +
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  xlab("") + ylab("False Discovery Rate") +ggtitle("10% DA Taxa (Rare Taxa)")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5))+ coord_flip()

Power_rare_avg <- Power_rare_results %>% group_by(Method) %>%
  summarise(Power_rare_avg = mean(Power_rare))  %>% arrange(Power_rare_avg)
Power_rare_results$Method <- factor(Power_rare_results$Method,
                                    levels=Power_rare_avg$Method)
Power_rare_plot <- ggplot(Power_rare_results, aes(x=Method, y=Power_rare)) + geom_violin(scale="count") +
  geom_boxplot(width=0.05) + ylim(0, 1) +
  xlab("") + ylab("Power(Sensitivity)") +ggtitle("10% DA Taxa (Rare Taxa)")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5))+ coord_flip()


FDR_abundant_avg <- FDR_abundant_results %>% group_by(Method) %>% 
  summarise(FDR_abundant_avg = mean(FDR_abundant)) %>% arrange(desc(FDR_abundant_avg))
FDR_abundant_results$Method <- factor(FDR_abundant_results$Method,
                                  levels=FDR_abundant_avg$Method)
FDR_abundant_plot <- ggplot(FDR_abundant_results, aes(x=Method, y=FDR_abundant)) + geom_violin(scale="count") +
  geom_boxplot(width=0.02) + ylim(0, 0.5) +
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  xlab("") + ylab("False Discovery Rate") +ggtitle("10% DA Taxa (Common Taxa)")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5))+ coord_flip()


Power_abundant_avg <- Power_abundant_results %>% group_by(Method) %>%
  summarise(Power_abundant_avg = mean(Power_abundant))  %>% arrange(Power_abundant_avg)
Power_abundant_results$Method <- factor(Power_abundant_results$Method,
                                    levels=Power_abundant_avg$Method)
Power_abundant_plot <- ggplot(Power_abundant_results, aes(x=Method, y=Power_abundant)) + geom_violin(scale="count") +
  geom_boxplot(width=0.05) + ylim(0, 1) +
  xlab("") + ylab("Power(Sensitivity)") +ggtitle("10% DA Taxa (Common Taxa)")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5))+ coord_flip()

library(patchwork)

(FDR_rare_plot + Power_rare_plot) / (FDR_abundant_plot + Power_abundant_plot)


auroc_combined <- cbind(aldex2_performance$AUROC,
                        ancom_performance$AUROC,
                        ancombc_performance$AUROC,
                        polda_performance$AUROC) |> as.data.frame()
colnames(auroc_combined) <- c("Aldex2", 'ANCOM', "ANCOMBC", "POLDA")
winner <- apply(auroc_combined, 1, which.max)

auroc_diff_df <- cbind(polda_performance$AUROC - aldex2_performance$AUROC,
                       polda_performance$AUROC - ancombc_performance$AUROC,
                       polda_performance$AUROC - ancom_performance$AUROC) |> as.data.frame()

colnames(auroc_diff_df) <- c("POLDA-ALDEx2", "POLDA-ANCOMBC", "POLDA-ANCOM")
auroc_diff_df_long <- auroc_diff_df %>% pivot_longer(cols=c("POLDA-ALDEx2", "POLDA-ANCOMBC", "POLDA-ANCOM"),
                                          names_to='Pair',
                                          values_to='Difference')



combined_performance <- rbind(aldex2_performance, ancom_performance, ancombc_performance, polda_performance)

library(cowplot)

AUROC_comparison_plot1 <- ggplot(combined_performance, aes(x=Method, y=AUROC)) + geom_violin()+
  geom_boxplot(width=0.02) + xlab("") + ylab("AUROC")+theme_bw()+coord_flip()

AUROC_comparison_plot2 <- ggplot(auroc_diff_df_long, aes(x=Pair, y=Difference)) + geom_violin()+
  geom_boxplot(width=0.02) +geom_hline(yintercept=0, linetype="dashed", color="red") +
  xlab("") + ylab("AUROC Difference")+theme_bw()+coord_flip()

plot_grid(AUROC_comparison_plot1, AUROC_comparison_plot2, ncol=1)

combined_performance_subset <- combined_performance %>% filter(Method != "ANCOM")

FDR_comparison_plot <- ggplot(combined_performance_subset, aes(x=Method, y=FDR)) + geom_violin() +
  geom_boxplot(width=0.02) + ylim(0, 0.4) +
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  xlab("") + ylab("False Discovery Rate") +
  theme_bw() + coord_flip()

power_comparison_plot <- ggplot(combined_performance_subset, aes(x=Method, y=Power)) + geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("") + ylab("Power") +
  theme_bw() + coord_flip()
library(cowplot)
plot_grid(FDR_comparison_plot, power_comparison_plot, ncol=1)

# POLDA_taxa_count$Falseref <- POLDA_taxa_count$Falseref > 0
# ggplot(POLDA_taxa_count, aes(x=RefTaxa, y=DiffTaxa, color=Falseref)) +
#   geom_jitter(size=1.5, alpha=0.7, width=0.1) + theme_bw() +
#   scale_color_manual(values=c("gray1", "brown2")) +
#   xlab("Number of Reference Taxa") + ylab("Number of DA Taxa")+
#   theme(legend.position = "None")
