rm(list=ls())
# performance_summary
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(pROC))

simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation'



IDs <- seq(1, 100)
modes <- c("mix")
foldchanges <- c(2.5, 3, 3.5, 4)
samplesizes <- c(100)
DAprops <- c(10)
seqdepths <- c(1e5)

create_df <- function(ID, modes, foldchange, samplesize, DAprop, seqdepth, method="PTDA"){
  performance_df <- expand.grid(ID, modes, foldchange, samplesize, DAprop, seqdepth)
  colnames(performance_df) <- c("IDs", "Modes", "FoldChanges", "SampleSizes", "DAProps", 'Seqdepths')
  performance_df$AUROC <- 0
  performance_df$FDR <- 0
  performance_df$Power <- 0
  performance_df$Time <- 0
  performance_df$Method <- method
  return(performance_df)
}

PTDA_folder <- file.path(simulation_folder, "PTDA")
PTDA_performance_noscale_df <- create_df(IDs, modes, foldchanges, samplesizes, DAprops, seqdepths, method="PTDA_noscale")
PTDA_performance_scale_df <- create_df(IDs, modes, foldchanges, samplesizes, DAprops, seqdepths, method="PTDA_scale")

  
for (rowID in 1:nrow(PTDA_performance_noscale_df)){
  ID <- PTDA_performance_noscale_df$IDs[rowID]
  mode <- PTDA_performance_noscale_df$Modes[rowID]
  foldchange <- PTDA_performance_noscale_df$FoldChanges[rowID]
  samplesize <- PTDA_performance_noscale_df$SampleSizes[rowID]
  DAprop <- PTDA_performance_noscale_df$DAProps[rowID]
  seqdepth <- PTDA_performance_noscale_df$Seqdepths[rowID]
  
  selected_mode <- PTDA_performance_noscale_df$Modes[rowID]
  input_filename <- sprintf("ptda_noboot_n%d_p%d_d%d_f%.1f_s%d.rds", 
                            samplesize, DAprop, seqdepth, foldchange, ID)
  model_summary <- readRDS(file.path(PTDA_folder, mode, sprintf("seed_%d", ID), input_filename))
  model_performance <- model_summary$performance
  PTDA_performance_noscale_df$Power[rowID] <- model_performance$Power
  PTDA_performance_noscale_df$AUROC[rowID] <- model_performance$AUROC
  PTDA_performance_noscale_df$FDR[rowID] <- model_performance$FDR
  PTDA_performance_noscale_df$Time[rowID] <- model_summary$duration[3]
}


for (rowID in 1:nrow(PTDA_performance_scale_df)){
  ID <- PTDA_performance_scale_df$IDs[rowID]
  mode <- PTDA_performance_scale_df$Modes[rowID]
  foldchange <- PTDA_performance_scale_df$FoldChanges[rowID]
  samplesize <- PTDA_performance_scale_df$SampleSizes[rowID]
  DAprop <- PTDA_performance_scale_df$DAProps[rowID]
  seqdepth <- PTDA_performance_scale_df$Seqdepths[rowID]
  
  selected_mode <- PTDA_performance_scale_df$Modes[rowID]
  input_filename <- sprintf("ptda_boot_n%d_p%d_d%d_f%.1f_s%d.rds", 
                            samplesize, DAprop, seqdepth, foldchange, ID)
  model_summary <- readRDS(file.path(PTDA_folder, mode, sprintf("seed_%d", ID), input_filename))
  model_performance <- model_summary$performance
  PTDA_performance_scale_df$Power[rowID] <- model_performance$Power
  PTDA_performance_scale_df$AUROC[rowID] <- model_performance$AUROC
  PTDA_performance_scale_df$FDR[rowID] <- model_performance$FDR
  PTDA_performance_scale_df$Time[rowID] <- model_summary$duration[3]
}


DACOMP_folder <- file.path(simulation_folder, "DACOMP")
DACOMP_performance_df <- create_df(IDs, modes, foldchanges, samplesizes, DAprops, seqdepths, method="DACOMP")

for (rowID in 1:nrow(DACOMP_performance_df)){
  ID <- DACOMP_performance_df$IDs[rowID]
  mode <- DACOMP_performance_df$Modes[rowID]
  foldchange <- DACOMP_performance_df$FoldChanges[rowID]
  samplesize <- DACOMP_performance_df$SampleSizes[rowID]
  DAprop <-DACOMP_performance_df$DAProps[rowID]
  seqdepth <- DACOMP_performance_df$Seqdepths[rowID]
  
  selected_mode <- DACOMP_performance_df$Modes[rowID]
  input_filename <- sprintf("dacomp_n%d_p%d_d%d_f%.1f_s%d.rds", 
                            samplesize, DAprop, seqdepth, foldchange, ID)
  
  model_summary<- readRDS(file.path(DACOMP_folder, mode, sprintf("seed_%d", ID), input_filename))
  model_performance <- model_summary$performance_dacomp
  DACOMP_performance_df$Power[rowID] <- model_performance$Power
  DACOMP_performance_df$AUROC[rowID] <- model_performance$AUROC
  DACOMP_performance_df$FDR[rowID] <- model_performance$FDR
  DACOMP_performance_df$Time[rowID] <- model_summary$duration[3]
}



ANCOMBC_folder <- file.path(simulation_folder, "ANCOMBC")
ANCOMBC_performance_df <- create_df(IDs, modes, foldchanges, samplesizes, DAprops, seqdepths, method="ANCOMBC")

for (rowID in 1:nrow(ANCOMBC_performance_df)){
  ID <- ANCOMBC_performance_df$IDs[rowID]
  mode <- ANCOMBC_performance_df$Modes[rowID]
  foldchange <- ANCOMBC_performance_df$FoldChanges[rowID]
  samplesize <- ANCOMBC_performance_df$SampleSizes[rowID]
  DAprop <-ANCOMBC_performance_df$DAProps[rowID]
  seqdepth <- ANCOMBC_performance_df$Seqdepths[rowID]
  
  selected_mode <- ANCOMBC_performance_df$Modes[rowID]
  input_filename <- sprintf("ANCOMBC_n%d_p%d_d%d_f%.1f_s%d.rds", 
                            samplesize, DAprop, seqdepth, foldchange, ID)
  model_summary <- readRDS(file.path(ANCOMBC_folder, mode, sprintf("seed_%d", ID), input_filename))
  model_performance <- model_summary$performance_ancombc
  ANCOMBC_performance_df$Power[rowID] <- model_performance$Power
  ANCOMBC_performance_df$AUROC[rowID] <- model_performance$AUROC
  ANCOMBC_performance_df$FDR[rowID] <- model_performance$FDR
  ANCOMBC_performance_df$Time[rowID] <- model_summary$duration[3]
}


# ANCOM_folder <- file.path(simulation_folder, "ANCOM")
# ANCOM_performance_df <- create_df(IDs, modes, foldchanges, samplesizes, DAprops, seqdepths, method="ANCOM")
# 
# 
# for (rowID in 1:nrow(ANCOM_performance_df)){
#   ID <- ANCOM_performance_df$IDs[rowID]
#   mode <- ANCOM_performance_df$Modes[rowID]
#   foldchange <- ANCOM_performance_df$FoldChanges[rowID]
#   samplesize <- ANCOM_performance_df$SampleSizes[rowID]
#   DAprop <-ANCOM_performance_df$DAProps[rowID]
#   seqdepth <- ANCOM_performance_df$Seqdepths[rowID]
#   selected_mode <- ANCOM_performance_df$Modes[rowID]
#   
#   input_filename <- sprintf("ANCOM_n%d_p%d_d%d_f%.1f_s%d.rds", 
#                             samplesize, DAprop, seqdepth, foldchange, ID)
#   model_performance <- readRDS(file.path(ANCOM_folder, mode, sprintf("seed_%d", ID), input_filename))
#   ANCOM_performance_df$Power[rowID] <- model_performance$Power
#   # ANCOM_performance_df$AUROC[rowID] <- model_performance$AUROC
#   ANCOM_performance_df$FDR[rowID] <- model_performance$FDR
# }


aldex_folder <- file.path(simulation_folder, "Aldex2")
aldex_performance_df <- create_df(IDs, modes, foldchanges, samplesizes, DAprops, seqdepths, method="ALDEx2")


for (rowID in 1:nrow(aldex_performance_df)){
  ID <- aldex_performance_df$IDs[rowID]
  mode <- aldex_performance_df$Modes[rowID]
  foldchange <- aldex_performance_df$FoldChanges[rowID]
  samplesize <- aldex_performance_df$SampleSizes[rowID]
  DAprop <-aldex_performance_df$DAProps[rowID]
  seqdepth <- aldex_performance_df$Seqdepths[rowID]
  selected_mode <- aldex_performance_df$Modes[rowID]
  input_filename <- sprintf("Aldex_n%d_p%d_d%d_f%.1f_s%d.rds", 
                            samplesize, DAprop, seqdepth, foldchange, ID)
  model_summary <- readRDS(file.path(aldex_folder, 
                                         mode, sprintf("seed_%d", ID), input_filename))
  model_performance <- model_summary$performance_aldex
  
  aldex_performance_df$Power[rowID] <- model_performance$Power
  aldex_performance_df$AUROC[rowID] <- model_performance$AUROC
  aldex_performance_df$FDR[rowID] <- model_performance$FDR
  aldex_performance_df$Time[rowID] <- model_summary$duration[3]
}



locom_folder <- file.path(simulation_folder, "LOCOM")
locom_performance_df <- create_df(IDs, modes, foldchanges, samplesizes, DAprops, seqdepths, method="LOCOM")


for (rowID in 1:nrow(locom_performance_df)){
  ID <- locom_performance_df$IDs[rowID]
  mode <- locom_performance_df$Modes[rowID]
  foldchange <- locom_performance_df$FoldChanges[rowID]
  samplesize <- locom_performance_df$SampleSizes[rowID]
  DAprop <-locom_performance_df$DAProps[rowID]
  seqdepth <- locom_performance_df$Seqdepths[rowID]
  selected_mode <- locom_performance_df$Modes[rowID]
  
  input_filename <- sprintf("LOCOM_n%d_p%d_d%d_f%.1f_s%d.rds", 
                            samplesize, DAprop, seqdepth, foldchange, ID)
  model_summary <- readRDS(file.path(locom_folder,
                                         mode, sprintf("seed_%d", ID), input_filename))
  model_performance <- model_summary$performance_locom
  locom_performance_df$Power[rowID] <- model_performance$Power
  locom_performance_df$AUROC[rowID] <- model_performance$AUROC
  locom_performance_df$FDR[rowID] <- model_performance$FDR
  locom_performance_df$Time[rowID] <- model_summary$duration[3]
}




all_results <- rbind(PTDA_performance_noscale_df,
                     PTDA_performance_scale_df,
                     aldex_performance_df,
                     ANCOMBC_performance_df,
                     DACOMP_performance_df,
                     locom_performance_df)

all_results$FoldChanges <- sprintf("Fold Change = %.1f", all_results$FoldChanges)

power_comparison <- ggplot(all_results, aes(x=Method, y=Power)) + geom_violin() + 
  geom_boxplot(width=0.1) + facet_grid(cols=vars(FoldChanges))+
  ylim(0, 1) +  ylab("Power")+theme_bw()+
  theme(text = element_text(size=10), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1))

FDR_comparison <- ggplot(all_results, aes(x=Method, y=FDR)) + geom_violin() + 
  geom_boxplot(width=0.1) + facet_grid(cols=vars(FoldChanges))+
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  ylim(0, 0.3) + ylab("FDR")+theme_bw()+
  theme(text = element_text(size=10), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1) )

library(cowplot)
plot_grid(FDR_comparison, power_comparison, ncol=1)

avg_times <- all_results %>% group_by(Method) %>%
  summarise(avgtime = mean(Time))

