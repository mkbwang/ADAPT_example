
rm(list=ls())
folder <- 'real_data/ECC/WGS'
comparison_plaque <- read.csv(file.path(folder, 
                                         "method_comparison_WGS_plaque.csv"))
colnames(comparison_plaque) <- c("Name", "ADAPT", "ALDEx2", "Maaslin2", "ZicoSeq",
                                             "ANCOM", "ANCOMBC", "LinDA", "metagenomeSeq", "DACOMP")
comparison_plaque[is.na(comparison_plaque)] <- 0
comparison_plaque <- comparison_plaque %>% dplyr::select(Name, ADAPT, ALDEx2, Maaslin2, metagenomeSeq, DACOMP,
                                                  ZicoSeq, ANCOM, ANCOMBC, LinDA)


library(reshape2)
comparison_WGS_long <- melt(comparison_plaque, id.vars="Name", variable.name="Method",
                           value.name="Direction")
comparison_WGS_long$Direction <- factor(comparison_WGS_long$Direction)
comparison_WGS_long$Name <- factor(comparison_WGS_long$Name, 
                                  levels=rev(comparison_plaque$Name))



library(ggplot2)
comparison_plot <- ggplot(comparison_WGS_long, aes(Method, Name, fill = Direction))+
  geom_tile(color="gray") + 
  scale_fill_manual(values=c("#077DE8", "white", '#F04520')) +
  xlab("Method") + ylab("Taxon") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none")






adapt_plaque_result <- read.csv(file.path(folder,
                                           "ADAPT_WGS_plaque_result.csv"))
adapt_plaque_result$neglog10pval <- -log10(adapt_plaque_result$pval)
adapt_plaque_result$log2effect <- adapt_plaque_result$effect * log2(exp(1))
union_taxa <- comparison_plaque$Taxon
adapt_taxa <- comparison_plaque$Taxon[comparison_plaque$ADAPT]
adapt_plaque_result <- adapt_plaque_result %>%
  rename(Taxon=Gene)


DAstatus <- as.integer(adapt_plaque_result$Taxon %in% union_taxa) +
  as.integer(adapt_plaque_result$Taxon %in% adapt_taxa)
DAtypes <- c("NonDA", "DA in competitors", "DA in ADAPT")
adapt_plaque_result$Type <-DAtypes[DAstatus +1]

adapt_DA_plaque_result <- adapt_plaque_result %>% 
  filter(Type=="DA in ADAPT")
adapt_nonDA_plaque_result <- adapt_plaque_result %>% 
  filter(Type!="DA in ADAPT")

logpval_cutoff <- (max(adapt_nonDA_plaque_result$neglog10pval) + 
                  min(adapt_DA_plaque_result$neglog10pval))/2

adapt_plaque_result$adapt_label <- NA
adapt_plaque_result$adapt_label[DAstatus==2] <- adapt_plaque_result$Taxon[DAstatus==2]
adapt_plaque_result$others_label <- NA
adapt_plaque_result$others_label[DAstatus == 1] <- adapt_plaque_result$Taxon[DAstatus == 1]

# volcano plot
library(ggplot2)
# library(ggrepel)
volcano_plaque <- ggplot(adapt_plaque_result, aes(x=log2effect, y=neglog10pval)) +
  geom_point(size=1.8, alpha=0.8, aes(color=Type)) + xlab("Log2FoldChange") + ylab("-log10pval") + theme_bw() + 
  theme(text = element_text(size=13), legend.position="bottom") + scale_color_manual(values=c("#ff0066", "#cc6600", "#666699")) +
  geom_hline(yintercept=logpval_cutoff, linetype="dashed", color = "blue") +
  geom_vline(xintercept=0, linetype="dashed", color = "blue")+
  scale_color_manual(values=c("#ff0066", "#cc6600", "#666699"))





# combine the plots
library(cowplot)
plot_grid(volcano_plaque, comparison_plot, ncol=2)


# complete dataframe
details_WGS_plaque <- adapt_plaque_result[, c("Taxon", "prevalence", 'log2effect', 'pval', 'adjusted_pval')] %>% 
  right_join(comparison_plaque)

write.csv(details_WGS_plaque, file.path(folder, "details_WGS_plaque_result.csv"),
          row.names = F)



