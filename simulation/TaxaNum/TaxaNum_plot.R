rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/TaxaNum"

settings_df <- expand.grid(nTaxa=c(100, 200, 500, 1000),
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
              sd_duration=sd(Duration, na.rm=T),
              Duration=mean(Duration, na.rm=T))
  
  # results_summary <- results_summary %>% filter(!Method %in% c("ADAPT_boot", "ADAPT_noboot", "LOCOM"))
  
  results_summary$nTaxa<- as.character(settings_df$nTaxa[choice])
  results_summary$Direction <- settings_df$direction[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)


all_summaries_df$nTaxa <- factor(all_summaries_df$nTaxa, 
                                  levels=c("100", "200", "500", "1000"))
all_summaries_df$Direction <- factor(all_summaries_df$Direction, 
                                     levels=c("Balanced Change", "Unbalanced Change"))
all_summaries_df$isADAPT <- all_summaries_df$Method == "ADAPT"



manual_color <- c("#666666", "#0066ff")
# Balanced 
subset_summaries_df_1 <- all_summaries_df %>% filter(Direction == "Balanced Change")
FDR_plot_1 <- ggplot(subset_summaries_df_1, aes(x=nTaxa, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.15))+
  scale_color_manual(values=manual_color)+
  xlab("Total Number of Taxa") + ylab("False Discovery Rate") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank())

FDR_plot_1

Power_plot_1 <- ggplot(subset_summaries_df_1, aes(x=nTaxa, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1))+
  scale_color_manual(values=manual_color)+
  xlab("Total Number of Taxa") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank())

Power_plot_1


# Unbalanced
subset_summaries_df_2 <- all_summaries_df %>% filter(Direction == "Unbalanced Change")
FDR_plot_2 <- ggplot(subset_summaries_df_2, aes(x=nTaxa, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.15))+
  scale_color_manual(values=manual_color)+
  xlab("Total Number of Taxa") + ylab("False Discovery Rate") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())

FDR_plot_2

Power_plot_2 <- ggplot(subset_summaries_df_2, aes(x=nTaxa, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1))+
  scale_color_manual(values=manual_color)+
  xlab("Total Number of Taxa") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())

Power_plot_2


write.csv(all_summaries_df, file.path(folder, "TaxaNum_summary.csv"),
          row.names=F)

time_df <- all_summaries_df %>% select(Method, sd_duration, Duration, nTaxa)
time_df$durations <- sprintf("%.4e(%.4e)", time_df$Duration, time_df$sd_duration)
time_df <- time_df %>% select(-c("Duration", "sd_duration"))
time_df$nTaxa <- as.character(time_df$nTaxa)
library(tidyr)
time_df_wide <- pivot_wider(time_df[seq(1, 36), ], id_cols="Method", names_from="nTaxa", values_from="durations")

write.csv(time_df_wide, file.path(folder, "TaxaNum_time.csv"),
          row.names=F)

# library(scales)
# duration_df <- all_summaries_df %>% group_by(Method, nTaxa) %>% summarise(Duration=mean(Duration))
# duration_df$isADAPT <- duration_df$Method == "ADAPT"
# time_plot <- ggplot(duration_df, aes(x=nTaxa, y=Duration, color=isADAPT, group=Method))+
#   geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
#   scale_y_continuous(trans='log10', 
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#   annotation_logticks(sides="l")+
#   scale_color_manual(values=manual_color)+
#   xlab("Total Number of Taxa") + ylab("Time (second)") + theme_bw()+
#   theme(text=element_text(size=14), legend.position = "None",
#         axis.title.x=element_blank(), axis.title.y=element_blank())
# 
# time_plot
