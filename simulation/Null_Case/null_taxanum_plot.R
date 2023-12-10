rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)


folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/Null_Case"

settings_df <- expand.grid(nSample=c(50, 100, 200),
                           nTaxa=c(100, 200, 500, 1000),
                           cf_main_corr=c(0, 0.3, 0.6),
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
    summarise(FPR=mean(FPR, na.rm=T),
              FWER=mean(FWER, na.rm=T),
              Duration=mean(Duration, na.rm=T))
  
  results_summary <- results_summary %>% filter(Method != "ADAPT_noboot") %>%
    mutate(Method=replace(Method, Method == "ADAPT_boot", "ADAPT"))
  
  results_summary$nSample <- as.character(settings_df$nSample[choice])
  results_summary$nTaxa <- settings_df$nTaxa[choice]
  results_summary$Corr <- settings_df$cf_main_corr[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)
subset_summaries <- all_summaries_df %>% filter(nSample == "100") %>% filter(!is.na(FPR))
subset_summaries$FPRLabel <- ""
subset_summaries$FWERLabel <- ""
subset_summaries$FPRLabel[subset_summaries$FPR > 0.1 & subset_summaries$nTaxa == "1000"] <- 
  subset_summaries$Method[subset_summaries$FPR > 0.1 & subset_summaries$nTaxa == "1000"]
subset_summaries$FPRLabel[subset_summaries$nTaxa == "1000" & subset_summaries$Method=="ADAPT"] <- "ADAPT"

subset_summaries$FWERLabel[subset_summaries$FWER > 0.1 & subset_summaries$nTaxa == "1000"] <- 
  subset_summaries$Method[subset_summaries$FWER > 0.1 & subset_summaries$nTaxa == "1000"]
subset_summaries$FWERLabel[subset_summaries$nTaxa == "1000" & subset_summaries$Method=="ADAPT"] <- "ADAPT"

subset_summaries$nTaxa <- factor(subset_summaries$nTaxa, 
                                      levels=c("100", "200", "500", "1000"))
subset_summaries$Corr <- sprintf("Corr(X, Z) = %.1f", subset_summaries$Corr)
subset_summaries$isADAPT <- subset_summaries$Method == "ADAPT"



manual_color <- c("#666666", "#d742f5")

FPR_plot <- ggplot(subset_summaries, aes(x=nTaxa, y=FPR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7, linetype="dashed") + 
  geom_hline(yintercept=0.05, linetype="dashed", color="red") +
  scale_y_continuous(limits=c(0, 0.25))+
  facet_grid(cols=vars(Corr)) + scale_color_manual(values=manual_color)+
  xlab("Taxa Number") + ylab("Type I Error") + theme_bw() + 
  geom_text_repel(aes(label = FPRLabel), nudge_x = 0.2,
                  na.rm = TRUE,  size=3.5)+
  theme(text=element_text(size=14), legend.position = "None")


FPR_plot

FWER_plot <- ggplot(subset_summaries, aes(x=nTaxa, y=FWER, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7, linetype="dashed") + 
  geom_hline(yintercept=0.05, linetype="dashed", color="red") +
  scale_y_continuous(limits=c(0, 1))+
  facet_grid(cols=vars(Corr)) + scale_color_manual(values=manual_color)+
  xlab("Taxa Number") + ylab("Familywise Error Rate") + theme_bw() + 
  geom_text_repel(aes(label = FWERLabel), nudge_x = 0.2,
                  na.rm = TRUE,  size=3.5)+
  theme(text=element_text(size=14), legend.position = "None")

FWER_plot

combined_plot <- plot_grid(FPR_plot, FWER_plot, nrow=2)

combined_plot


