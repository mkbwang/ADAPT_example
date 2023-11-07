
rm(list=ls())
folder <- 'real_data/ECC/16S'
comparison_24month <- read.csv(file.path(folder, 
                                         "method_comparison_24month.csv"),
                               row.names = 1)


adapt_24month_result <- read.csv(file.path(folder,
                                           "ADAPT_24month_result.csv"))
adapt_24month_result$neglog10pval <- -log10(adapt_24month_result$pval)
adapt_24month_result$log2effect <- adapt_24month_result$effect * log2(exp(1))
union_taxa <- comparison_24month$Taxon
adapt_taxa <- comparison_24month$Taxon[comparison_24month$ADAPT]
adapt_24month_result <- adapt_24month_result %>%
  rename(Taxon=Gene)


DAstatus <- as.integer(adapt_24month_result$Taxon %in% union_taxa) +
  as.integer(adapt_24month_result$Taxon %in% adapt_taxa)
DAtypes <- c("NonDA", "DA in competitors", "DA in ADAPT")
adapt_24month_result$Type <-DAtypes[DAstatus +1]

adapt_DA_24month_result <- adapt_24month_result %>% 
  filter(Type=="DA in ADAPT")
adapt_nonDA_24month_result <- adapt_24month_result %>% 
  filter(Type!="DA in ADAPT")

logpval_cutoff <- (max(adapt_nonDA_24month_result$neglog10pval) + 
                  min(adapt_DA_24month_result$neglog10pval))/2

adapt_24month_result$adapt_label <- NA
adapt_24month_result$adapt_label[DAstatus==2] <- adapt_24month_result$Taxon[DAstatus==2]
adapt_24month_result$others_label <- NA
adapt_24month_result$others_label[DAstatus == 1] <- adapt_24month_result$Taxon[DAstatus == 1]

# volcano plot
library(ggplot2)
# library(ggrepel)
volcano_24month <- ggplot(adapt_24month_result, aes(x=log2effect, y=neglog10pval)) +
  geom_point(size=1.8, alpha=0.8, aes(color=Type)) + xlab("Log2FoldChange") + ylab("-log10pval") + theme_bw() + 
  theme(text = element_text(size=13)) + scale_color_manual(values=c("#ff0066", "#cc6600", "#666699")) +
  geom_hline(yintercept=logpval_cutoff, linetype="dashed", color = "blue") +
  geom_vline(xintercept=0, linetype="dashed", color = "blue")+
  scale_color_manual(values=c("#ff0066", "#cc6600", "#666699"))



# heatmap plot of DAA result comparison
DA_full_df <- comparison_24month[, seq(1, 9)]
library(reshape2)
DA_long_df <- melt(DA_full_df, id.vars="Taxon", variable.name="Method",
                   value.name="isDA")
DA_long_df$isDA <- as.integer(DA_long_df$isDA)
DA_long_df$isDA[DA_long_df$Method == "ADAPT" & DA_long_df$isDA] <- 2
DA_long_df$isDA <- as.factor(DA_long_df$isDA)
DA_long_df$Taxon <- factor(DA_long_df$Taxon, levels=DA_long_df$Taxon[1:30])


comparison_plot <- ggplot(DA_long_df, aes(Taxon, Method, fill = isDA))+
  geom_tile(color="gray") + 
  scale_fill_manual(values=c("white", "black", '#3366ff')) +
  xlab("Taxon") + ylab("Method") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none")

# combine the plots
library(cowplot)
plot_grid(volcano_24month, comparison_plot, ncol=1)


# complete dataframe
details_16S_24month <- adapt_24month_result[, c("Taxon", "prevalence", 'log2effect', 'pval', 'adjusted_pval')] %>% 
  right_join(comparison_24month)

write.csv(details_16S_24month, file.path(folder, "details_24month_result.csv"),
          row.names = F)



