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

results_ADAPT <- all_results_df %>% filter(Method == "ADAPT") %>%
  filter(depth_fold != "Unbalanced Library Size (4-fold)")
results_ADAPT$depth_fold <- as.character(results_ADAPT$depth_fold)
results_ADAPT$depth_fold[which(results_ADAPT$depth_fold == "Unbalanced Library Size (10-fold)")] = "Unbalanced Library Size"

decision_counts_adapt <- results_ADAPT %>% group_by(nSample, depth_fold, num_FD) %>%
  summarise(Count=n())

decision_counts_adapt$num_FD <- as.character(decision_counts_adapt$num_FD)

ggplot(decision_counts_adapt, aes(x=num_FD, y=Count)) + geom_bar(stat="identity", color="black", fill="lightblue")+
  facet_grid(rows=vars(depth_fold), cols=vars(nSample)) +
  xlab("Number of False Discoveries") + ylab("Frequency") + theme_bw() + coord_flip()+ 
  theme(text=element_text(size=12),
        strip.text = element_text(size = 10))






# all_results_df$Method[all_results_df$Method == "Maaslin2"] <- "MaAsLin2" 

all_summary <- all_results_df %>% group_by(Method, nSample, depth_fold) %>%
  summarise(FPR_1=mean(FPR_1, na.rm=T), FPR_5 = mean(FPR_5, na.rm=T), FPR_10=mean(FPR_10, na.rm=T),
            num_FD = mean(num_FD, na.rm=T))


manual_color <- c("#666666", "#0066ff")



# remove the cases corresponding to 4-fold sequencing depth difference
subset_results_df <- all_results_df %>% filter(depth_fold != "Unbalanced Library Size (4-fold)")
subset_results_df$depth_fold <- as.character(subset_results_df$depth_fold)
subset_results_df$depth_fold[subset_results_df$depth_fold == "Unbalanced Library Size (10-fold)"] = "Unbalanced Library Size"

FPR_5_plot <- ggplot(subset_results_df, aes(x=Method, y=FPR_5, color=isADAPT, group=Method)) +
  geom_boxplot() +
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.6), breaks=seq(0, 0.6, 0.1))+
  facet_grid(rows=vars(depth_fold), cols=vars(nSample)) + scale_color_manual(values=manual_color)+
  xlab("") + ylab("") + theme_bw() + 
  theme(text=element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title=element_blank(),
        legend.position = "None",
        strip.text = element_text(size = 10))

FPR_5_plot




write.csv(all_summary, file.path(folder, "Null_performances.csv"),
          row.names=F)




