
rm(list=ls())
folder <- 'real_data/ECC/WGS'
## load the DAA result of all nine methods
comparison_plaque <- read.csv(file.path(folder, 
                                         "method_comparison_WGS_plaque.csv"))
## load the details of ADAPT output
adapt_plaque_result <- read.csv(file.path(folder,
                                          "ADAPT_WGS_plaque_result.csv"))
adapt_plaque_result$Taxa <- gsub("[.]" , "", adapt_plaque_result$Taxa)
phy_metag_plaque <- readRDS(file.path(folder, "phy_metag_plaque.rds"))
comparison_plaque$Taxa <- trimws(comparison_plaque$Taxa)

# First generate the supplementary table
comparison_plaque <- comparison_plaque %>% dplyr::select(Taxa, ADAPT,
                                                  ALDEx2, MaAsLin2, metagenomeSeq, DACOMP,
                                                  ZicoSeq,ANCOM,ANCOMBC, LinDA)
detailed_metag_plaque <- comparison_plaque %>% 
  full_join(adapt_plaque_result, by="Taxa")
write.csv(detailed_metag_plaque, file.path(folder, "detailed_WGS_plaque_result.csv"),
          row.names=F)


# visualize the DAA result of different methods
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



# volcano plot of the ADAPT result
adapt_plaque_result$neglog10pval <- -log10(adapt_plaque_result$pval)
adapt_plaque_result$log10effect <- adapt_plaque_result$effect * log10(exp(1))
union_taxa <- comparison_plaque$Taxa
adapt_taxa <- comparison_plaque$Taxa[comparison_plaque$ADAPT!=0]



DAstatus <- as.integer(adapt_plaque_result$Taxa %in% union_taxa) +
  as.integer(adapt_plaque_result$Taxa %in% adapt_taxa)
DAtypes <- c("NonDA", "DA in competitors", "DA in ADAPT")
adapt_plaque_result$Type <-DAtypes[DAstatus +1]


adapt_plaque_result$adapt_label <- NA
adapt_plaque_result$adapt_label[DAstatus==2] <- adapt_plaque_result$Taxa[DAstatus==2]
adapt_plaque_result$others_label <- NA
adapt_plaque_result$others_label[DAstatus == 1] <- adapt_plaque_result$Taxa[DAstatus == 1]


library(ggplot2)
# library(ggrepel)
volcano_plaque <- ggplot(adapt_plaque_result, aes(x=log10effect, y=neglog10pval)) +
  geom_point(alpha=0.8, aes(color=Type)) + 
  xlab("Log10 Fold Change") + ylab("-Log10 p-value") + theme_bw() + 
  theme(legend.position="none", axis.title=element_blank(), axis.text=element_text(size=10)) + scale_color_manual(values=c("#ff0066", "#cc6600", "#666699")) +
  geom_vline(xintercept=0, linetype="dashed", color = "blue")+
  scale_x_continuous(breaks=seq(-9, 10, 1))+
  scale_y_continuous(breaks=seq(0, 8, 1))


