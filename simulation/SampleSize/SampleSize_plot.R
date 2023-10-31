rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/SampleSize"

settings_df <- expand.grid(Direction=c("balanced", "unbalanced"),
                           Nsample=c(50, 100, 200))


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
  
  results_summary$Nsample <- as.character(settings_df$Nsample[choice])
  results_summary$Direction <- settings_df$Direction[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)
all_summaries_df$Nsample <- factor(all_summaries_df$Nsample, 
                                  levels=c("50", "100", "200"))
all_summaries_df$Direction <- factor(all_summaries_df$Direction, 
                                     levels=c("balanced", "unbalanced"))


FDR_plot <- ggplot(all_summaries_df, aes(x=Nsample, y=FDR, color=Method, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=1, alpha=0.7) + 
  geom_hline(yintercept=0.05, linetype="dashed", color="red") +
  facet_grid(cols=vars(Direction)) + scale_color_brewer(palette="Dark2")+
  xlab("Sample Size") + ylab("False Discovery Rate") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "bottom")

Power_plot <- ggplot(all_summaries_df, aes(x=Nsample, y=Power, color=Method, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=1, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1))+
  facet_grid(cols=vars(Direction)) + scale_color_brewer(palette="Dark2")+
  xlab("Sample Size") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "bottom")


combined_plot <- plot_grid(FDR_plot, Power_plot, nrow=1)




