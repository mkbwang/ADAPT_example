
rm(list=ls())
folder <- 'real_data/ECC/WGS'
comparison_plaque <- read.csv(file.path(folder, 
                                         "method_comparison_WGS_plaque.csv"))
comparison_plaque[is.na(comparison_plaque)] <- 0


library(reshape2)
comparison_WGS_long <- melt(comparison_plaque, id.vars="Taxa", variable.name="Method",
                           value.name="Direction")
comparison_WGS_long$Direction <- factor(comparison_WGS_long$Direction)
comparison_WGS_long$Taxa <- factor(comparison_WGS_long$Taxa, 
                                  levels=rev(comparison_plaque$Taxa))



library(ggplot2)
comparison_plot <- ggplot(comparison_WGS_long, aes(Method, Taxa, fill = Direction))+
  geom_tile(color="gray") + 
  scale_fill_manual(values=c("#077DE8", "white", '#F04520')) +
  xlab("Method") + ylab("Taxon") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none",
        axis.title = element_blank())






adapt_plaque_result <- read.csv(file.path(folder,
                                           "ADAPT_WGS_plaque_result.csv"))
adapt_plaque_result$Taxa <- gsub("[.]", " ", adapt_plaque_result$Taxa)
adapt_plaque_result$neglog10pval <- -log10(adapt_plaque_result$pval)
adapt_plaque_result$log2effect <- adapt_plaque_result$effect * log2(exp(1))
union_taxa <- comparison_plaque$Taxa
adapt_taxa <- comparison_plaque$Taxa[comparison_plaque$ADAPT!=0]



DAstatus <- as.integer(adapt_plaque_result$Taxa %in% union_taxa) +
  as.integer(adapt_plaque_result$Taxa %in% adapt_taxa)
DAtypes <- c("NonDA", "DA in competitors", "DA in ADAPT")
adapt_plaque_result$Type <-DAtypes[DAstatus +1]

# adapt_DA_plaque_result <- adapt_plaque_result %>% 
#   filter(Type=="DA in ADAPT")
# adapt_nonDA_plaque_result <- adapt_plaque_result %>% 
#   filter(Type!="DA in ADAPT")
# 
# logpval_cutoff <- (max(adapt_nonDA_plaque_result$neglog10pval) + 
#                   min(adapt_DA_plaque_result$neglog10pval))/2

adapt_plaque_result$adapt_label <- NA
adapt_plaque_result$adapt_label[DAstatus==2] <- adapt_plaque_result$Taxa[DAstatus==2]
adapt_plaque_result$others_label <- NA
adapt_plaque_result$others_label[DAstatus == 1] <- adapt_plaque_result$Taxa[DAstatus == 1]

# volcano plot
library(ggplot2)
# library(ggrepel)
volcano_plaque <- ggplot(adapt_plaque_result, aes(x=log2effect, y=neglog10pval)) +
  geom_point(alpha=0.8, aes(color=Type)) + 
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") + theme_bw() + 
  theme(legend.position="none") + scale_color_manual(values=c("#ff0066", "#cc6600", "#666699")) +
  geom_vline(xintercept=0, linetype="dashed", color = "blue")+
  scale_x_continuous(breaks=seq(-18, 20, 2))+
  scale_y_continuous(breaks=seq(0, 8, 1))







# complete dataframe
details_WGS_plaque <- adapt_plaque_result[, c("Taxa", "prevalence", 'log2effect', 'pval', 'adjusted_pval')] %>% 
  right_join(comparison_plaque)

write.csv(details_WGS_plaque, file.path(folder, "details_WGS_plaque_result.csv"),
          row.names = F)



