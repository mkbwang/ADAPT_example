rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/Dimension"

settings_df <- expand.grid(nSample=c(50, 100, 200),
                           nTaxa=c(100, 200, 500, 1000),
                           direction=c("Balanced Change", "Unbalanced Change"),
                           stringsAsFactors = FALSE)


all_summaries <- list()

files <- list.files(file.path(folder, 'experiments'))

choices <- seq(1, nrow(settings_df))

for (choice in choices){
  
  name_filter <- grepl(sprintf("experiment_%d", choice), 
                       files)
  subset_files <- files[name_filter]
  all_results <- lapply(subset_files, 
                        function(fname) readRDS(file.path(folder, "experiments", fname)))
  
  all_results_df <- do.call(rbind, all_results)
  results_summary <- all_results_df %>% group_by(Method) %>%
    summarise(FDR=mean(FDR, na.rm=T),
              Power=mean(Power, na.rm=T),
              Duration=mean(Duration, na.rm=T))
  
  results_summary <- results_summary %>% filter(Method != "ADAPT_noboot") %>%
    mutate(Method=replace(Method, Method == "ADAPT_boot", "ADAPT"))
  
  results_summary$SampleSize<- as.character(settings_df$nSample[choice])
  results_summary$TaxaNum <- as.character(settings_df$nTaxa[choice])
  results_summary$Direction <- settings_df$direction[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)
subset_summaries <- all_summaries_df %>% filter(SampleSize == "100")

subset_summaries$FDRLabel <- ""
subset_summaries$PowerLabel <- ""
subset_summaries$FDRLabel[subset_summaries$FDR > 0.055 & subset_summaries$TaxaNum == "100"] <- 
  subset_summaries$Method[subset_summaries$FDR > 0.055 & subset_summaries$TaxaNum == "100"]
subset_summaries$FDRLabel[subset_summaries$TaxaNum == "100" & subset_summaries$Method=="ADAPT"] <- "ADAPT"

subset_summaries$PowerLabel[subset_summaries$TaxaNum == "100"] <- 
  subset_summaries$Method[subset_summaries$TaxaNum == "100"]

subset_summaries$TaxaNum <- factor(subset_summaries$TaxaNum, 
                                  levels=c("100", "200", "500", "1000"))
subset_summaries$Direction <- factor(subset_summaries$Direction, 
                                     levels=c("Balanced Change", "Unbalanced Change"))
subset_summaries$isADAPT <- subset_summaries$Method == "ADAPT"



manual_color <- c("#666666", "#d742f5")

FDR_plot <- ggplot(subset_summaries, aes(x=TaxaNum, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7, linetype="dashed") + 
  geom_hline(yintercept=0.05, linetype="dashed", color="red") +
  scale_y_continuous(limits=c(0, 0.1))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Taxa Number") + ylab("False Discovery Rate") + theme_bw() + 
  geom_text_repel(aes(label = FDRLabel), nudge_x =- 0.2,
                  na.rm = TRUE,  size=3.5)+
  theme(text=element_text(size=14), legend.position = "None")


FDR_plot

Power_plot <- ggplot(subset_summaries, aes(x=TaxaNum, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7, linetype="dashed") + 
  scale_y_continuous(limits=c(0, 1))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Taxa Number") + ylab("Power") + theme_bw() + 
  geom_text_repel(aes(label = PowerLabel),
                  nudge_x = -0.2,
                  na.rm = TRUE,  size=3.5)+
  theme(text=element_text(size=14), legend.position = "None")


Power_plot

combined_plot <- plot_grid(FDR_plot, Power_plot, nrow=1)

combined_plot

