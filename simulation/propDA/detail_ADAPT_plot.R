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

files <- list.files(file.path(folder, 'details_ADAPT'))



choices <- seq(1, nrow(settings_df))
for (choice in choices){
  
  name_filter <- grepl(sprintf("SparseDOSSA_%d_", choice), # SparseDOSSA
                       files)
  subset_files <- files[name_filter]
  all_results <- lapply(subset_files, 
                        function(fname) readRDS(file.path(folder, "details_ADAPT", fname)))
  
  all_results_df <- do.call(rbind, lapply(all_results, as.data.frame))
  
  # results_summary <- all_results_df %>% group_by(Method) %>%
  #   summarise(
  #     sd_FDR = sd(FDR, na.rm=T),
  #     FDR=mean(FDR, na.rm=T),
  #     sd_Power = sd(Power, na.rm=T),
  #     Power = mean(Power, na.rm=T),
  #     Duration=mean(Duration, na.rm=T))
  
  all_results_df$propDA<- as.character(settings_df$propDA[choice]*100)
  all_results_df$Direction <- settings_df$direction[choice]
  all_summaries[[choice]] <- all_results_df
  
}



all_summaries_df <- do.call(rbind, all_summaries)
all_summaries_df$propDA <- factor(all_summaries_df$propDA, 
                                  levels=c("5", "10", "20", "30"))
all_summaries_df$Direction <- factor(all_summaries_df$Direction, 
                                     levels=c("Balanced Change", "Unbalanced Change"))
all_summaries_df$IoU <- all_summaries_df$num_intersection / all_summaries_df$num_union
all_summaries_df$refprop <- sprintf("1/%d", 
                                    round(all_summaries_df$num_taxa / all_summaries_df$num_reftaxa)) 
all_summaries_df$contamination <- all_summaries_df$ref_error / all_summaries_df$num_reftaxa * 100

# contamination proportion
manual_color <- c("#3438eb", "#c934eb")
ggplot(all_summaries_df, aes(x=propDA, y=contamination, color=refprop)) + geom_jitter(alpha=0.5) + 
  xlab("Proportion of DA Taxa among All Taxa (%)") + ylab("Proportion of DA Taxa among Reference Taxa (%)")+
  facet_wrap(~Direction) + scale_color_manual(values=manual_color)+
  theme_bw() + theme(text=element_text(size=12))


# Overlap between detected taxa and alternative detected taxa
ggplot(all_summaries_df, aes(x=propDA, y=IoU)) + geom_violin() + geom_boxplot(width=0.1)+
  xlab("Proportion of DA Taxa among All Taxa (%)") + ylab("Jaccard Similarity Coefficient")+
  facet_wrap(~Direction) + 
  theme_bw() + theme(text=element_text(size=12))


all_summaries_df$FDRdiff <- all_summaries_df$FDR - all_summaries_df$altFDR
all_summaries_df$Powerdiff <- all_summaries_df$power - all_summaries_df$altPower
# Boxplot of FDR difference and power difference 
ggplot(all_summaries_df, aes(x=propDA, y=FDRdiff)) + geom_violin() + geom_boxplot(width=0.1)+
  xlab("Proportion of DA Taxa among All Taxa (%)") + ylab("Difference of FDR")+
  facet_wrap(~Direction) + 
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  theme_bw() + theme(text=element_text(size=12))


ggplot(all_summaries_df, aes(x=propDA, y=Powerdiff)) + geom_violin() + geom_boxplot(width=0.1)+
  xlab("Proportion of DA Taxa among All Taxa (%)") + ylab("Difference of Power")+
  facet_wrap(~Direction) + 
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  theme_bw() + theme(text=element_text(size=12))

all_summaries_df$isADAPT <- all_summaries_df$Method == "ADAPT"



