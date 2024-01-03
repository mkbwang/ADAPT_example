rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/propDA"

settings_df <- expand.grid(propDA=c(0.05, 0.1, 0.2, 0.3),
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
  
  results_summary$propDA<- as.character(settings_df$propDA[choice]*100)
  results_summary$Direction <- settings_df$direction[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)
all_summaries_df$propDA <- factor(all_summaries_df$propDA, 
                                  levels=c("5", "10", "20", "30"))
all_summaries_df$Direction <- factor(all_summaries_df$Direction, 
                                     levels=c("Balanced Change", "Unbalanced Change"))
all_summaries_df$isADAPT <- all_summaries_df$Method == "ADAPT"



manual_color <- c("#666666", "#0066ff")

FDR_plot <- ggplot(all_summaries_df, aes(x=propDA, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.35))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Proportion of DA Taxa (%)") + ylab("False Discovery Rate") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None")


FDR_plot

Power_plot <- ggplot(all_summaries_df, aes(x=propDA, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Proportion of DA Taxa (%)") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None")


Power_plot

write.csv(all_summaries_df, file.path(folder, "propDA_summary.csv"),
          row.names=F)

