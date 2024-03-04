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

  
  results_summary$SampleSize<- as.character(settings_df$nSample[choice])
  results_summary$Direction <- settings_df$direction[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)
all_summaries_df <- all_summaries_df %>% filter(SampleSize != "30")
all_summaries_df$SampleSize <- factor(all_summaries_df$SampleSize, 
                                  levels=c( "50", "80", "100", "150", "200"))
all_summaries_df$Direction <- factor(all_summaries_df$Direction, 
                                     levels=c("Balanced Change", "Unbalanced Change"))
all_summaries_df$isADAPT <- all_summaries_df$Method == "ADAPT"



# balanced, FDR
manual_color <- c("#666666", "#0066ff")
subset_summaries_df_1 <- all_summaries_df %>% filter(Direction == "Balanced Change")
FDR_plot_1 <- ggplot(subset_summaries_df_1, aes(x=SampleSize, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.2, 0.05))+
  scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("False Discovery Rate") + theme_bw() +
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank())


FDR_plot_1


# unbalanced, FDR
subset_summaries_df_2 <- all_summaries_df %>% filter(Direction == "Unbalanced Change")
FDR_plot_2 <- ggplot(subset_summaries_df_2, aes(x=SampleSize, y=FDR, color=isADAPT, group=Method)) +
  geom_point(size=1.4, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dotted", color="red") +
  scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.2, 0.05))+
  scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("False Discovery Rate") + theme_bw() +
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())


FDR_plot_2



# balanced, Power
Power_plot_1 <- ggplot(subset_summaries_df_1, aes(x=SampleSize, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1))+
  scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank())

Power_plot_1

# unbalanced, Power
Power_plot_2 <- ggplot(subset_summaries_df_2, aes(x=SampleSize, y=Power, color=isADAPT, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1))+
  scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())

Power_plot_2


write.csv(all_summaries_df, 
          file.path(folder, "SampleSize_summary.csv"),
          row.names=F)


# the duration of all methods with respect to the sample size (not used)
duration_df <- all_summaries_df %>% group_by(Method, SampleSize) %>%
  summarise(Duration=mean(Duration))
duration_df$isADAPT <- duration_df$Method == "ADAPT"
library(scales)
time_plot <- ggplot(duration_df, aes(x=SampleSize, y=Duration, color=isADAPT, group=Method))+
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=0.8, alpha=0.7) + 
  scale_y_continuous(trans='log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="l")+
  scale_color_manual(values=manual_color)+
  xlab("Sample Size") + ylab("Time (second)") + theme_bw()+
  theme(text=element_text(size=14), legend.position = "None",
        axis.title.x=element_blank(), axis.title.y=element_blank())

time_plot
