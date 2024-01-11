rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/Null_Case"

settings_df <- expand.grid(nSample=c("Sample Size = 50", "Sample Size = 100"),
                           nTaxa=500,
                           depth_fold=c("Balanced Library Size", "Unbalanced Library Size (4-fold)",
                                        "Unbalanced Library Size (10-fold)"))


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
  results_df$depth_fold <- settings_df$depth_fold[choice]
  all_results[[choice]] <- results_df
}

all_results_df <- do.call(rbind, all_results)
# subset_results <- all_results_df %>% filter(Method != "DACOMP")
all_results_df$isADAPT <- all_results_df$Method == "ADAPT"
all_results_df$nSample <- factor(all_results_df$nSample,
                                 levels=c("Sample Size = 50", "Sample Size = 100"))
all_results_df$depth_fold <- factor(all_results_df$depth_fold,
                                  levels=c("Balanced Library Size", "Unbalanced Library Size (4-fold)",
                                           "Unbalanced Library Size (10-fold)"))

all_results_df$Method[all_results_df$Method == "Maaslin2"] <- "MaAsLin2" 

all_summary <- all_results_df %>% group_by(Method, nSample, depth_fold) %>%
  summarise(FPR_1=mean(FPR_1, na.rm=T), FPR_5 = mean(FPR_5, na.rm=T), FPR_10=mean(FPR_10, na.rm=T),
            FWE_1 = mean(FWE_1, na.rm=T))


manual_color <- c("#666666", "#0066ff")

# FPR_1_noadj_plot <- ggplot(subset_results_noadj, aes(x=Method, y=FPR_1, color=isADAPT, group=Method)) +
#   geom_boxplot() +
#   geom_hline(yintercept=0.01, linetype="dotted", color="red") +
#   scale_y_continuous(limits=c(0, 0.15))+
#   facet_grid(cols=vars(balance_depth)) + scale_color_manual(values=manual_color)+
#   xlab("Method") + ylab("False Positive Rates at level of 0.01") + theme_bw() + 
#   theme(text=element_text(size=14),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         legend.position = "None")


# FPR_1_noadj_plot


FPR_1_plot <- ggplot(all_results_df, aes(x=Method, y=FPR_1, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.01, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.4))+
  facet_grid(cols=vars(depth_fold), rows=vars(nSample)) + scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.01") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")

FPR_1_plot


FPR_5_plot <- ggplot(all_results_df, aes(x=Method, y=FPR_5, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.6), breaks=seq(0, 0.6, 0.1))+
  facet_grid(cols=vars(depth_fold), rows=vars(nSample)) + scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rate") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")

FPR_5_plot



FPR_10_plot <- ggplot(all_results_df, aes(x=Method, y=FPR_10, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.10, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.8))+
  facet_grid(rows=vars(depth_fold), vars=vars(nSample)) + scale_color_manual(values=manual_color)+
  xlab("Method") + ylab("False Positive Rates at level of 0.1") + theme_bw() + 
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "None")

FPR_10_plot




write.csv(all_summary, file.path(folder, "Null_performances.csv"),
          row.names=F)




