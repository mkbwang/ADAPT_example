rm(list=ls())
# performance_summary
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)

data_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/data"
old_POLDA_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/POLDA_nobiascorrect"
new_POLDA_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/POLDA"
locom_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/LOCOM"
metagenomeseq_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/MetagenomeSeq"

polda_old_performance <- readRDS(file.path(old_POLDA_folder, 
                                           "polda_old_performance.rds"))

polda_new_performance <- readRDS(file.path(new_POLDA_folder, 
                                           "polda_new_performance.rds"))

locom_performance <- readRDS(file.path(locom_folder,
                                       "locom_performance.rds"))

metagenomeseq_performance <- readRDS(file.path(metagenomeseq_folder,
                                               "metagenomeseq_performance.rds"))

all_performances <- rbind(polda_old_performance,
                          polda_new_performance,
                          locom_performance,
                          metagenomeseq_performance)


performances_100 <- all_performances %>% filter(sample_size==100)
performances_100$Proportion <- factor(as.character(performances_100$Proportion),
                                      levels=c("5", "10", "20"))

performances_50 <- all_performances %>% filter(sample_size == 50)
performances_50$Proportion <- factor(as.character(performances_50$Proportion),
                                      levels=c("5", "10", "20"))


power_comparison_100 <- ggplot(performances_100, aes(x=Proportion, y=Power, color=Method)) + 
  geom_boxplot(width=0.6) + facet_grid(rows=vars(Direction), cols=vars(Type))+
  scale_color_manual(values=c("brown", "gray1", "blue", "darkgreen"))+ ggtitle("Sample Size 100")+
  ylim(0, 1.01) + xlab("DA Proportion (%)") + ylab("Power")+theme_bw()+
  theme(text = element_text(size=13), plot.title = element_text(hjust = 0.5))

FDR_comparison_100 <- ggplot(performances_100, aes(x=Proportion, y=FDR, color=Method)) + 
  geom_boxplot(width=0.6) + facet_grid(rows=vars(Direction), cols=vars(Type))+
  scale_color_manual(values=c("brown", "gray1", "blue", "darkgreen"))+
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+ ggtitle("Sample Size 100")+
  ylim(0, 0.45) + xlab("DA Proportion (%)") + ylab("FDR")+theme_bw()+
  theme(text = element_text(size=14), plot.title = element_text(hjust = 0.5))

AUROC_comparison_100 <- ggplot(performances_100, aes(x=Proportion, y=AUROC, color=Method)) + 
  geom_boxplot(width=0.6) + facet_grid(rows=vars(Direction), cols=vars(Type))+
  scale_color_manual(values=c("brown", "gray1", "blue", "darkgreen"))+
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+ ggtitle("Sample Size 100")+
  ylim(0.8, 1) + xlab("DA Proportion (%)") + ylab("AUROC")+theme_bw()+
  theme(text = element_text(size=14), plot.title = element_text(hjust = 0.5))


power_comparison_50 <- ggplot(performances_50, aes(x=Proportion, y=Power, color=Method)) + 
  geom_boxplot(width=0.6) + facet_grid(rows=vars(Direction), cols=vars(Type))+
  scale_color_manual(values=c("brown", "gray1", "blue", "darkgreen"))+ ggtitle("Sample Size 50")+
  ylim(0, 1.01) + xlab("DA Proportion (%)") + ylab("Power")+theme_bw()+
  theme(text = element_text(size=14), plot.title = element_text(hjust = 0.5))

FDR_comparison_50 <- ggplot(performances_50, aes(x=Proportion, y=FDR, color=Method)) + 
  geom_boxplot(width=0.6) + facet_grid(rows=vars(Direction), cols=vars(Type))+
  scale_color_manual(values=c("brown", "gray1", "blue", "darkgreen"))+ ggtitle("Sample Size 50")+
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  ylim(0, 0.45) + xlab("DA Proportion (%)") + ylab("FDR")+theme_bw()+
  theme(text = element_text(size=14), plot.title = element_text(hjust = 0.5))

AUROC_comparison_50 <- ggplot(performances_50, aes(x=Proportion, y=AUROC, color=Method)) + 
  geom_boxplot(width=0.6) + facet_grid(rows=vars(Direction), cols=vars(Type))+
  scale_color_manual(values=c("brown", "gray1", "blue", "darkgreen"))+
  ggtitle("Sample Size 50")+
  ylim(0.8, 1) + xlab("DA Proportion (%)") + ylab("AUROC")+theme_bw()+
  theme(text = element_text(size=14), plot.title = element_text(hjust = 0.5))


