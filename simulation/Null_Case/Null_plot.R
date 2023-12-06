rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/Null_Case"

settings_df <- expand.grid(nSample=c(50, 100),
                           balance_depth = c("Unbalanced Library Size", "Balanced Library Size"),
                           adj_depth = c("NoLibAdj", "LibAdj"))


all_results <- list()

files <- list.files(file.path(folder, 'experiments'))

choices <- seq(1, nrow(settings_df))
for (choice in choices){
  
  name_filter <- grepl(sprintf("experiment_%d_", choice), 
                       files)
  subset_files <- files[name_filter]
  results <- lapply(subset_files, 
                        function(fname) readRDS(file.path(folder, "experiments", fname)))
  
  results_df <- do.call(rbind, results)
  
  # results_summary <- results_summary %>% filter(!Method %in% c("ADAPT_boot", "ADAPT_noboot", "LOCOM"))
  results_df$nSample <- settings_df$nSample[choice]
  results_df$balance_depth<- settings_df$balance_depth[choice]
  results_df$adj_depth <- settings_df$adj_depth[choice]
  all_results[[choice]] <- results_df
  
}

all_results_df <- do.call(rbind, all_results)
subset_results_noadj <- all_results_df %>% filter(nSample==100 & adj_depth == "NoLibAdj")
subset_results_noadj$isADAPT <- subset_results_noadj$Method == "ADAPT"
subset_results_noadj$balance_depth <- factor(subset_results_noadj$balance_depth,
                                             levels=c("Balanced Library Size", "Unbalanced Library Size"))
subset_noadj_summary <- subset_results_noadj %>% group_by(Method, balance_depth) %>%
  summarise(FPR_1=median(FPR_1, na.rm=T), FPR_5 = median(FPR_5, na.rm=T), FPR_10=median(FPR_10, na.rm=T))


manual_color <- c("#666666", "#0066ff")

FPR_1_noadj_plot <- ggplot(subset_results_noadj, aes(x=Method, y=FPR_1, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.01, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.15))+
  facet_grid(cols=vars(balance_depth)) + scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.01") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")


FPR_1_noadj_plot

FPR_5_noadj_plot <- ggplot(subset_results_noadj, aes(x=Method, y=FPR_5, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.3))+
  facet_grid(cols=vars(balance_depth)) + scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.05") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")


FPR_5_noadj_plot

FPR_10_noadj_plot <- ggplot(subset_results_noadj, aes(x=Method, y=FPR_10, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.1, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.4))+
  facet_grid(cols=vars(balance_depth)) + scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.1") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")


FPR_10_noadj_plot


subset_results_adj <- all_results_df %>% filter(nSample==100 & adj_depth == "LibAdj" & balance_depth=="Unbalanced Library Size") %>% 
  filter(Method %in% c("ADAPT", "ALDEx2", "ANCOMBC", "LinDA", "Maaslin2", "ZicoSeq"))
subset_results_adj$isADAPT <- subset_results_adj$Method == "ADAPT"
subset_results_adj$balance_depth <- factor(subset_results_adj$balance_depth,
                                             levels=c("Balanced Library Size", "Unbalanced Library Size"))
subset_adj_summary <- subset_results_adj %>% group_by(Method, balance_depth) %>%
  summarise(FPR_1=median(FPR_1, na.rm=T), FPR_5 = median(FPR_5, na.rm=T), FPR_10=median(FPR_10, na.rm=T))



FPR_1_adj_plot <- ggplot(subset_results_adj, aes(x=Method, y=FPR_1, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.01, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.1))+
  scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.01") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")


FPR_1_adj_plot


FPR_5_adj_plot <- ggplot(subset_results_adj, aes(x=Method, y=FPR_5, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.25))+
  scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.05") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")


FPR_5_adj_plot


FPR_10_adj_plot <- ggplot(subset_results_adj, aes(x=Method, y=FPR_10, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.1, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.35))+
  scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.1") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")


FPR_10_adj_plot

write.csv(subset_noadj_summary, file.path(folder, "Null_LibSize_Noadjust_performances.csv"),
          row.names=F)

write.csv(subset_adj_summary, file.path(folder, "Null_LibSize_adjust_performances.csv"),
          row.names=F)

