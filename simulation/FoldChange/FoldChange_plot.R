rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/FoldChange"

settings_df <- expand.grid(FoldChange=c(3, 3.5, 4, 5, 6),
  direction=c("Balanced Change", "Unbalanced Change"),
                           stringsAsFactors = FALSE)


all_summaries <- list()

files <- list.files(file.path(folder, 'experiments'))

choices <- seq(1, nrow(settings_df))
for (choice in choices){
  
  name_filter <- grepl(sprintf("experiment_%d_", choice), 
                       files)
  subset_files <- files[name_filter]
  all_results <- lapply(subset_files, 
                        function(fname) readRDS(file.path(folder, "experiments", fname)))
  
  all_results_df <- do.call(rbind, all_results)
  results_summary <- all_results_df %>% group_by(Method) %>%
    summarise(FDR=mean(FDR, na.rm=T),
              Power=mean(Power, na.rm=T),
              Duration=mean(Duration, na.rm=T))
  
  # results_summary <- results_summary %>% filter(!Method %in% c("ADAPT_boot", "ADAPT_noboot", "LOCOM"))
  
  results_summary$FoldChange<- as.character(settings_df$FoldChange[choice])
  results_summary$Direction <- settings_df$direction[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)

# all_summaries_df$FDRLabel <- ""
# all_summaries_df$PowerLabel <- ""
# all_summaries_df$FDRLabel[all_summaries_df$FDR > 0.1] <- 
#   all_summaries_df$Method[all_summaries_df$FDR > 0.1]
# all_summaries_df$FDRLabel[all_summaries_df$FoldChange=="3" & all_summaries_df$Method=="ADAPT"] <- "ADAPT"

# all_summaries_df$PowerLabel[all_summaries_df$FoldChange=="6"] <- 
#   all_summaries_df$Method[all_summaries_df$FoldChange=="6"]
# 
all_summaries_df$FoldChange <- factor(all_summaries_df$FoldChange,
                                  levels=c("3", "3.5", "4", "5", "6"))
all_summaries_df$Direction <- factor(all_summaries_df$Direction,
                                     levels=c("Balanced Change", "Unbalanced Change"))
all_summaries_df$isADAPT <- all_summaries_df$Method == "ADAPT"



manual_color <- c("#666666", "#0066ff")

FDR_plot <- ggplot(all_summaries_df, aes(x=FoldChange, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.15))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Average Fold Change") + ylab("False Discovery Rate") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None")


FDR_plot

Power_plot <- ggplot(all_summaries_df, aes(x=FoldChange, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Average Fold Change") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None")


Power_plot


write.csv(all_summaries_df, file.path(folder,  "FoldChange_Summary.csv"),
          row.names=F)


