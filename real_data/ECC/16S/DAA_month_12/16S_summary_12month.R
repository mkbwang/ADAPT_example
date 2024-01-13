library(dplyr)
rm(list=ls())
folder <- 'real_data/ECC/16S'
comparison_12month <- read.csv(file.path(folder, "DAA_month_12",
                                         "method_comparison_12month.csv"))

Taxon_12 <- mapply(function(gname, sname) paste(gname, sname, collapse=" "),
                comparison_12month$Genus, comparison_12month$Species)

# heatmap plot of DAA result comparison
comparison_12month$Name <- Taxon_12
comparison_12month$Name <- factor(comparison_12month$Name,
                                     levels=unique(comparison_12month$Name))

## there are multiple ASVs with the same species name
comparison_12month_deduplicate <- comparison_12month %>% select(Name, ADAPT,
                                                                ALDEx2, MaAsLin2, metagenomeSeq, DACOMP,
                                                                ZicoSeq,ANCOM,ANCOMBC, LinDA) %>%
  group_by(Name) %>% summarise(across(everything(), sum))

comparison_12month_deduplicate[comparison_12month_deduplicate == -2] <- -1
comparison_12month_deduplicate[comparison_12month_deduplicate == 2] <- 1

library(reshape2)
comparison_12_long <- melt(comparison_12month_deduplicate, id.vars="Name", variable.name="Method",
                   value.name="Direction")
comparison_12_long$Direction <- factor(comparison_12_long$Direction)
comparison_12_long$Name <- factor(comparison_12_long$Name, 
                                  levels=rev(comparison_12month_deduplicate$Name))

library(ggplot2)
comparison_plot <- ggplot(comparison_12_long, aes(Method, Name, fill = Direction))+
  geom_tile(color="gray") + 
  scale_fill_manual(values=c("#077DE8", "white", '#F04520')) +
  xlab("Method") + ylab("Taxon") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none",
        axis.title = element_blank())




# volcano plot
adapt_12month_result <- read.csv(file.path(folder, "DAA_month_12",
                                           "ADAPT_12month_result.csv"))
adapt_12month_result$neglog10pval <- -log10(adapt_12month_result$pval)
adapt_12month_result$log2effect <- adapt_12month_result$effect * log2(exp(1))
union_taxa <- comparison_12month$Taxa
adapt_taxa <- comparison_12month$Taxa[comparison_12month$ADAPT != 0]


DAstatus <- as.integer(adapt_12month_result$Taxa %in% union_taxa) +
  as.integer(adapt_12month_result$Taxa %in% adapt_taxa)
DAtypes <- c("NonDA", "DA in competitors", "DA in ADAPT")
adapt_12month_result$Type <-DAtypes[DAstatus +1]

# adapt_DA_12month_result <- adapt_12month_result %>% 
#   filter(Type=="DA in ADAPT")
# adapt_nonDA_12month_result <- adapt_12month_result %>% 
#   filter(Type!="DA in ADAPT")
# 
# logpval_cutoff <- (max(adapt_nonDA_12month_result$neglog10pval) + 
#                   min(adapt_DA_12month_result$neglog10pval))/2

adapt_12month_result$adapt_label <- NA
adapt_12month_result$adapt_label[DAstatus==2] <- adapt_12month_result$Taxa[DAstatus==2]
adapt_12month_result$others_label <- NA
adapt_12month_result$others_label[DAstatus == 1] <- adapt_12month_result$Taxa[DAstatus == 1]

library(ggplot2)
# library(ggrepel)
volcano_12month <- ggplot(adapt_12month_result, aes(x=log2effect, y=neglog10pval)) +
  geom_point(alpha=0.8, aes(color=Type)) + 
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") + theme_bw() + 
  theme(legend.position="none") + scale_color_manual(values=c("#ff0066", "#cc6600", "#666699")) +
  geom_vline(xintercept=0, linetype="dashed", color = "blue")+
  scale_x_continuous(breaks=seq(-18, 20, 2))+
  scale_y_continuous(breaks=seq(0, 8, 1))







# complete dataframe
details_16S_12month <- adapt_12month_result[, c("Taxa", "prevalence", 'log2effect', 'pval', 'adjusted_pval')] %>% 
  right_join(comparison_12month, by="Taxa")

write.csv(details_16S_12month, file.path(folder, "details_12month_result.csv"),
          row.names = F)




