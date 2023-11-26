rm(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/FWER"

settings_df <- expand.grid(Depth_confound=c(TRUE, FALSE), 
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
    summarise(FWER=mean(FWE, na.rm=T))
  
  results_summary <- results_summary %>% filter(Method != "ADAPT_noboot") %>%
    mutate(Method=replace(Method, Method == "ADAPT_boot", "ADAPT"))
  
  results_summary$Nsample <- as.character(settings_df$Nsample[choice])
  results_summary$DepthConfound <- settings_df$Depth_confound[choice]
  all_summaries[[choice]] <- results_summary
  
}

all_summaries_df <- do.call(rbind, all_summaries)
all_summaries_df$Nsample <- factor(all_summaries_df$Nsample, 
                                  levels=c("50", "100", "200"))
all_summaries_df$DepthConfound <- factor(all_summaries_df$DepthConfound)

all_summaries_df_subset <- all_summaries_df %>% filter(DepthConfound == "FALSE")


FWER_plot <- ggplot(all_summaries_df_subset, aes(x=Nsample, y=FWER, color=Method, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=1, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 0.1))+
  geom_hline(yintercept=0.05, linetype="dashed", color="red") +
  scale_color_brewer(palette="Dark2")+
  xlab("Sample Size") + ylab("Familywise Error Rate") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "bottom")

Power_plot <- ggplot(all_summaries_df, aes(x=Nsample, y=Power, color=Method, group=Method)) +
  geom_point(size=1.2, alpha=0.8) + geom_line(linewidth=1, alpha=0.7) + 
  scale_y_continuous(limits=c(0, 1))+
  facet_grid(cols=vars(Direction)) + scale_color_brewer(palette="Dark2")+
  xlab("Sample Size") + ylab("Power") + theme_bw() + 
  theme(text=element_text(size=14), legend.position = "bottom")


combined_plot <- plot_grid(FDR_plot, Power_plot, nrow=1)




