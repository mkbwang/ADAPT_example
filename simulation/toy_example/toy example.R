# absolute abundance

aa1 <- data.frame(Taxon=sprintf("Taxon %d", seq(1, 7)),
                  Count=c(160, 100, 80, 40, 40, 10, 10))

aa2 <- data.frame(Taxon=sprintf("Taxon %d", seq(1, 7)),
                  Count=c(100, 100, 80, 0, 40, 10, 20))

library(ggplot2)
library(ggpattern)
library(RColorBrewer)
brewer.pal(n = 8, name = "Dark2")

ggplot(aa1, aes(x=factor(Taxon), y=Count, fill=Taxon)) + geom_bar(stat="identity") + 
  coord_flip() + xlab("") + 
  scale_y_continuous(name="Absolute Abundance", limits=c(0, 180), breaks=seq(0, 160, 20))+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none")


ggplot(aa2, aes(x=factor(Taxon), y=Count, fill=Taxon)) + geom_bar(stat="identity", color="000000") + 
  coord_flip() + xlab("") + ylab("Absolute Abundance") +
  scale_fill_brewer(palette = "Dark2") + scale_y_continuous(name="Absolute Abundance", limits=c(0, 180), breaks=seq(0, 160, 20))+
  theme_classic() + theme(legend.position = "none")


# sequenced metagenomic counts

ma1 <- data.frame(Taxon=sprintf("Taxon %d", seq(1, 7)),
                  Count=c(8, 5, 4, 2, 2, 0, 0))

ma2 <- data.frame(Taxon=sprintf("Taxon %d", seq(1, 7)),
                  Count=c(10, 10, 8, 0, 4, 1, 2))

ggplot(ma1, aes(x=factor(Taxon), y=Count, fill=Taxon)) + geom_bar(stat="identity", color="000000") + 
  coord_flip() + xlab("") + 
  scale_y_continuous(name="Sequencing Read Count", limits=c(0, 14), breaks=seq(0, 12, 2))+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none")




ggplot(ma2, aes(x=factor(Taxon), y=Count, fill=Taxon)) + geom_bar(stat="identity", color="000000") + 
  coord_flip() + xlab("") + 
  scale_y_continuous(name="Sequencing Read Count", limits=c(0, 14), breaks=seq(0, 12, 2))+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none")


# Relative Abundance


ra1 <- ma1
ra1$Count <- ra1$Count / sum(ra1$Count)
ggplot(ra1, aes(x=factor(Taxon), y=Count, fill=Taxon)) + geom_bar(stat="identity", color="000000") + 
  coord_flip() + xlab("") + 
  scale_y_continuous(name="Relative Abundance", limits=c(0, 0.4), breaks=seq(0, 8, 2)/21,
                     labels=c("0",  "2/21", "4/21",  "2/7",  "8/21"))+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none")

ra2 <- ma2
ra2$Count <- ra2$Count / sum(ra2$Count)
ggplot(ra2, aes(x=factor(Taxon), y=Count, fill=Taxon)) + geom_bar(stat="identity", color="000000") + 
  coord_flip() + xlab("") + 
  scale_y_continuous(name="Relative Abundance", limits=c(0, 0.4), breaks=seq(0, 8, 2)/21,
                     labels=c("0",  "2/21", "4/21",  "2/7",  "8/21"))+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none")



# Relative Abundance with censoring

ra1_copy <- ra1
ra1_copy$Censored <- "No"
ra1_copy$Censored[ra1_copy$Count == 0] <- "Yes"
ra1_copy$Count[ra1_copy$Count == 0] <- 1/sum(ma1$Count)
ggplot(ra1_copy, aes(x=factor(Taxon), y=Count, fill=Taxon, pattern=Censored)) +
  geom_bar_pattern(stat="identity", pattern_fill = "#222222",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
  coord_flip() + xlab("") + 
  scale_y_continuous(name="Relative Abundance", limits=c(0, 0.4), breaks=seq(0, 8, 2)/21,
                     labels=c("0",  "2/21", "4/21",  "2/7",  "8/21"))+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none")



ra2_copy <- ra2
ra2_copy$Censored <- "No"
ra2_copy$Censored[ra2_copy$Count == 0] <- "Yes"
ra2_copy$Count[ra2_copy$Count == 0] <- 1/sum(ma2$Count)
ggplot(ra2_copy, aes(x=factor(Taxon), y=Count, fill=Taxon, pattern=Censored)) +
  geom_bar_pattern(stat="identity", pattern_fill = "#222222",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
  coord_flip() + xlab("") + 
  scale_y_continuous(name="Relative Abundance", limits=c(0, 0.4), breaks=seq(0, 8, 2)/21,
                     labels=c("0",  "2/21", "4/21",  "2/7",  "8/21"))+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none")



# Scale by the reference set

cr1 <- ma1
refcounts1=sum(cr1$Count[c(2,3,5)])
cr1$Count <- cr1$Count / refcounts1
cr1$Censored <- "No"
cr1$Censored[cr1$Count == 0] <- "Yes"
cr1$Count[cr1$Count == 0] <- 1/refcounts1
cr1$Sample <- "Sample 1"
cr1$DA <- FALSE
cr1$DA[c(1,4,7)] <- TRUE 



cr2 <- ma2
refcounts2=sum(cr2$Count[c(2,3,5)])
cr2$Count <- cr2$Count / refcounts2
cr2$Censored <- "No"
cr2$Censored[cr2$Count == 0] <- "Yes"
cr2$Count[cr2$Count == 0] <- 1/refcounts2
cr2$Sample <- "Sample 2"
cr2$DA <- FALSE
cr2$DA[c(1, 4, 7)] <- TRUE

cr_combined <- rbind(cr1, cr2)
ggplot(cr_combined, aes(x=factor(Taxon), y=Count, fill=Sample, color=DA, pattern=Censored)) +
  geom_bar_pattern(stat="identity", pattern_fill = "#222222",
                   position=position_dodge(),
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
  scale_color_manual(values=c("#000000", "#BB0000"))+
  scale_fill_manual(values=c("#FFFFFF", "#999999"))+
  xlab("") +
  scale_y_continuous(name="Normalized Counts", limits=c(0, 0.8), breaks=seq(0, 8, 2)/11)+
  theme_classic() + theme(legend.position = "none", axis.title=element_blank(), axis.text=element_blank())

