rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/SampleSize"

settings_df <- expand.grid(nSample=c(30, 50, 80, 100, 150, 200),
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
  
  results_summary$SampleSize<- as.character(settings_df$nSample[choice])
  results_summary$Direction <- settings_df$direction[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)
all_summaries_df <- all_summaries_df %>% filter(SampleSize != "30")
# all_summaries_df$FDRLabel <- ""
# all_summaries_df$PowerLabel <- ""
# FDRlabel_mask <- all_summaries_df$FDR > 0.1 & (all_summaries_df$SampleSize == "50" | all_summaries_df$SampleSize == "200")
# all_summaries_df$FDRLabel[FDRlabel_mask] <- 
#   all_summaries_df$Method[FDRlabel_mask]
# all_summaries_df$FDRLabel[all_summaries_df$SampleSize == "50" & all_summaries_df$Method=="ADAPT"] <- "ADAPT"
# 
# all_summaries_df$PowerLabel[all_summaries_df$SampleSize=="200"] <- 
#   all_summaries_df$Method[all_summaries_df$SampleSize=="200"]

all_summaries_df$SampleSize <- factor(all_summaries_df$SampleSize, 
                                  levels=c( "50", "80", "100", "150", "200"))
all_summaries_df$Direction <- factor(all_summaries_df$Direction, 
                                     levels=c("Balanced Change", "Unbalanced Change"))
all_summaries_df$isADAPT <- all_summaries_df$Method == "ADAPT"




manual_color <- c("#666666", "#0066ff")

FDR_plot <- ggplot(all_summaries_df, aes(x=SampleSize, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.2))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("False Discovery Rate") + theme_bw() +
  theme(text=element_text(size=14), legend.position = "None")


FDR_plot

Power_plot <- ggplot(all_summaries_df, aes(x=SampleSize, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1))+
  facet_grid(cols=vars(Direction)) + scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None")


Power_plot

write.csv(all_summaries_df, 
          file.path(folder, "SampleSize_summary.csv"),
          row.names=F)

duration_df <- all_summaries_df %>% group_by(Method, SampleSize) %>%
  summarise(Duration=mean(Duration))
duration_df$isADAPT <- duration_df$Method == "ADAPT"
library(scales)
time_plot <- ggplot(duration_df, aes(x=SampleSize, y=Duration, color=isADAPT, group=Method))+
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(trans='log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("Time (second)") + theme_bw()+
  theme(text=element_text(size=14), legend.position = "None")

time_plot
