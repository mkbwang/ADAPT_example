library(dplyr)
rm(list=ls())
folder <- 'real_data/ECC/16S'
comparison_24month <- read.csv(file.path(folder, "DAA_month_24",
                                         "method_comparison_24month.csv"))

comparison_24month[is.na(comparison_24month)] <- "Unknown"

Taxon_24 <- mapply(function(gname, sname) paste(gname, sname, collapse=" "),
                   comparison_24month$Genus, comparison_24month$Species)
Taxon_24[30] <- "Unknown"

# heatmap plot of DAA result comparison
comparison_24month$Name <- Taxon_24
comparison_24month$Name <- factor(comparison_24month$Name,
                                  levels=unique(comparison_24month$Name))

## there are multiple ASVs with the same species name
comparison_24month_deduplicate <- comparison_24month %>% select(Name, ADAPT,
                                                                ALDEx2, MaAsLin2, metagenomeSeq, DACOMP,
                                                                ZicoSeq,ANCOM,ANCOMBC, LinDA) %>%
  group_by(Name) %>% summarise(across(everything(), sum))

comparison_24month_deduplicate[comparison_24month_deduplicate == -2] <- -1
comparison_24month_deduplicate[comparison_24month_deduplicate == 2] <- 1



library(reshape2)
comparison_24_long <- melt(comparison_24month_deduplicate, id.vars="Name", variable.name="Method",
                           value.name="Direction")
comparison_24_long$Direction <- factor(comparison_24_long$Direction)
comparison_24_long$Name <- factor(comparison_24_long$Name, 
                                  levels=rev(comparison_24month_deduplicate$Name))

library(ggplot2)
comparison_plot <- ggplot(comparison_24_long, aes(Method, Name, fill = Direction))+
  geom_tile(color="gray") + 
  scale_fill_manual(values=c("#077DE8", "white", '#F04520')) +
  xlab("Method") + ylab("Taxon") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none")



# volcano plot
adapt_24month_result <- read.csv(file.path(folder, "DAA_month_24",
                                           "ADAPT_24month_result.csv"))
adapt_24month_result$neglog10pval <- -log10(adapt_24month_result$pval)
adapt_24month_result$log2effect <- adapt_24month_result$effect * log2(exp(1))
union_taxa <- comparison_24month$Taxa
adapt_taxa <- comparison_24month$Taxa[comparison_24month$ADAPT != 0]



DAstatus <- as.integer(adapt_24month_result$Taxa %in% union_taxa) +
  as.integer(adapt_24month_result$Taxa %in% adapt_taxa)
DAtypes <- c("NonDA", "DA in competitors", "DA in ADAPT")
adapt_24month_result$Type <-DAtypes[DAstatus +1]

# adapt_DA_24month_result <- adapt_24month_result %>% 
#   filter(Type=="DA in ADAPT")
# adapt_nonDA_24month_result <- adapt_24month_result %>% 
#   filter(Type!="DA in ADAPT")
# 
# logpval_cutoff <- (max(adapt_nonDA_24month_result$neglog10pval) + 
#                   min(adapt_DA_24month_result$neglog10pval))/2

adapt_24month_result$adapt_label <- NA
adapt_24month_result$adapt_label[DAstatus==2] <- adapt_24month_result$Taxa[DAstatus==2]
adapt_24month_result$others_label <- NA
adapt_24month_result$others_label[DAstatus == 1] <- adapt_24month_result$Taxa[DAstatus == 1]

# volcano plot
library(ggplot2)
# library(ggrepel)
volcano_24month <- ggplot(adapt_24month_result, aes(x=log2effect, y=neglog10pval)) +
  geom_point(alpha=0.8, aes(color=Type)) + 
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") + theme_bw() + 
  theme(legend.position="none") + scale_color_manual(values=c("#ff0066", "#cc6600", "#666699")) +
  geom_vline(xintercept=0, linetype="dashed", color = "blue")+
  scale_x_continuous(breaks=seq(-18, 20, 2))+
  scale_y_continuous(breaks=seq(0, 8, 1))





# complete dataframe
details_16S_24month <- adapt_24month_result[, c("Taxa", "prevalence", 'log2effect', 'pval', 'adjusted_pval')] %>% 
  right_join(comparison_24month, by="Taxa")

write.csv(details_16S_24month, file.path(folder, "details_24month_result.csv"),
          row.names = F)


