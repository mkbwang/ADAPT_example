library(dplyr)
rm(list=ls())
folder <- 'real_data/ECC/16S'


comparison_saliva <- read.csv(file.path(folder, "method_comparison_16S_saliva_12month.csv"))
adapt_saliva_result <- read.csv(file.path(folder, "ADAPT_16S_saliva_12month_result.csv"))
phy_asv_saliva <- readRDS(file.path(folder, "phyloseq", "phyasv_visit12.rds"))
taxonomies <- data.frame(tax_table(phy_asv_saliva))
adapt_saliva_result <- cbind(adapt_saliva_result, 
                              taxonomies[adapt_saliva_result$Taxa, ])
all_DNAstrings <- as.character(phy_asv_saliva@refseq)
all_DNAstrings <- all_DNAstrings[adapt_saliva_result$Taxa]
adapt_saliva_result$DNAstring <- all_DNAstrings
prev_counts <- colSums(otu_table(phy_asv_saliva) > 0)

# First generate the supplementary table
detailed_adapt_saliva_result <- comparison_saliva[, seq(1, 10)] %>%
  right_join(adapt_saliva_result, by="Taxa")
detailed_adapt_saliva_result$prevalence <- prev_counts[detailed_adapt_saliva_result$Taxa]
sample_size <- nrow(sample_data(phy_asv_saliva))
prev_counts <- prev_counts[detailed_adapt_saliva_result$Taxa]
detailed_adapt_saliva_result$prevalence <- sprintf("%d/%d", prev_counts, sample_size)
detailed_adapt_saliva_result$effect <- detailed_adapt_saliva_result$effect * log10(exp(1)) # change natural log to log10
write.csv(detailed_adapt_saliva_result, file.path(folder, "detailed_16S_saliva_12month_result.csv"),
          row.names=F)


taxon_names <- mapply(function(gname, sname) paste(gname, sname, collapse=" "),
                comparison_saliva$Genus, comparison_saliva$Species)

# heatmap plot of DAA result comparison
comparison_saliva$Name <- taxon_names
comparison_saliva$Name <- factor(comparison_saliva$Name,
                                     levels=unique(comparison_saliva$Name))

## there are multiple ASVs with the same species name
comparison_saliva_deduplicate <- comparison_saliva %>% dplyr::select(Name, ADAPT,
                                                                ALDEx2, MaAsLin2, metagenomeSeq, DACOMP,
                                                                ZicoSeq,ANCOM,ANCOMBC, LinDA) %>%
  group_by(Name) %>% summarise(across(everything(), sum))

comparison_saliva_deduplicate[comparison_saliva_deduplicate == -2] <- -1
comparison_saliva_deduplicate[comparison_saliva_deduplicate == 2] <- 1

library(reshape2)
comparison_saliva_long <- melt(comparison_saliva_deduplicate, id.vars="Name", variable.name="Method",
                   value.name="Direction")
comparison_saliva_long$Direction <- factor(comparison_saliva_long$Direction)
comparison_saliva_long$Name <- factor(comparison_saliva_long$Name, 
                                  levels=rev(comparison_saliva_deduplicate$Name))

library(ggplot2)
comparison_plot <- ggplot(comparison_saliva_long, aes(Method, Name, fill = Direction))+
  geom_tile(color="gray") + 
  scale_fill_manual(values=c("#077DE8", "white", '#F04520')) +
  xlab("Method") + ylab("Taxon") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none",
        axis.title = element_blank())




# volcano plot
adapt_saliva_result$neglog10pval <- -log10(adapt_saliva_result$pval)
adapt_saliva_result$log10effect <- adapt_saliva_result$effect * log10(exp(1))
union_taxa <- comparison_saliva$Taxa
adapt_taxa <- comparison_saliva$Taxa[comparison_saliva$ADAPT != 0]


DAstatus <- as.integer(adapt_saliva_result$Taxa %in% union_taxa) +
  as.integer(adapt_saliva_result$Taxa %in% adapt_taxa)
DAtypes <- c("NonDA", "DA in competitors", "DA in ADAPT")
adapt_saliva_result$Type <-DAtypes[DAstatus +1]


adapt_saliva_result$adapt_label <- NA
adapt_saliva_result$adapt_label[DAstatus==2] <- adapt_saliva_result$Taxa[DAstatus==2]
adapt_saliva_result$others_label <- NA
adapt_saliva_result$others_label[DAstatus == 1] <- adapt_saliva_result$Taxa[DAstatus == 1]

library(ggplot2)
volcano_saliva <- ggplot(adapt_saliva_result, aes(x=log10effect, y=neglog10pval)) +
  geom_point(alpha=0.8, aes(color=Type)) + 
  xlab("Log10 Fold Change") + ylab("-Log10 p-value") + theme_bw() + 
  theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text=element_text(size=10)) + scale_color_manual(values=c("#ff0066", "#cc6600", "#666699")) +
  geom_vline(xintercept=0, linetype="dashed", color = "blue")+
  scale_x_continuous(breaks=seq(-4, 7, 1))+
  scale_y_continuous(breaks=seq(0, 8, 1))


