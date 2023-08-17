library(ggplot2)
simulation_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation"

cox_performance <- read.table(file.path(simulation_folder, "FPR", "permutation", "cox_FPR.txt"),
                              header=F, sep='\t')
colnames(cox_performance) <- c("Run", "FPR_5", "FPR_10", "FPR_20", "FPR_30")
cox_performance$Model <- "Cox"


glm_performance <- read.table(file.path(simulation_folder, "FPR", "permutation", "GLM_FPR.txt"),
                              header=F, sep='\t')
colnames(glm_performance) <- c("Run", "FPR_5", "FPR_10", "FPR_20", "FPR_30")
glm_performance$Model <- "GLM"


lm_performance <- read.table(file.path(simulation_folder, "FPR", "permutation", "lm_FPR.txt"),
                              header=F, sep='\t')
colnames(lm_performance) <- c("Run", "FPR_5", "FPR_10", "FPR_20", "FPR_30")
lm_performance$Model <- "Linear"


wilcox_performance <- read.table(file.path(simulation_folder, "FPR", "permutation", "wilcox_FPR.txt"),
                             header=F, sep='\t')
colnames(wilcox_performance) <- c("Run", "FPR_5", "FPR_10", "FPR_20", "FPR_30")
wilcox_performance$Model <- "Wilcox"


loglogistic_performance <- read.table(file.path(simulation_folder, "FPR", "permutation", "loglogistic_FPR.txt"),
                                 header=F, sep='\t')
colnames(loglogistic_performance) <- c("Run", "FPR_5", "FPR_10", "FPR_20", "FPR_30")
loglogistic_performance$Model <- "Loglogistic"



weibull_performance <- read.table(file.path(simulation_folder, "FPR", "permutation", "weibull_FPR.txt"),
                                      header=F, sep='\t')
colnames(weibull_performance) <- c("Run", "FPR_5", "FPR_10", "FPR_20", "FPR_30")
weibull_performance$Model <- "Weibull"


all_performances <- rbind(cox_performance, lm_performance, wilcox_performance, loglogistic_performance,
                          weibull_performance, glm_performance)
all_performances$Model <- factor(all_performances$Model, levels=c("Cox", "Loglogistic", "Weibull", "GLM", "Linear", "Wilcox"))


cutoff_5 <- ggplot(all_performances, aes(x=Model, y=FPR_5)) + geom_violin() + geom_boxplot(width=0.1) + 
  xlab("Method") + ylab("Type I Error") + ggtitle("Cutoff at 0.05") + 
  geom_hline(yintercept=0.05, linetype="dashed", color="red", linewidth=1)+ylim(0, 0.3)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=14))

cutoff_10 <- ggplot(all_performances, aes(x=Model, y=FPR_10)) + geom_violin() + geom_boxplot(width=0.1) + 
  xlab("Method") + ylab("Type I Error") + ggtitle("Cutoff at 0.1") + 
  geom_hline(yintercept=0.1, linetype="dashed", color="red", linewidth=1)+ylim(0, 0.3)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=14))


cutoff_20 <- ggplot(all_performances, aes(x=Model, y=FPR_20)) + geom_violin() +  geom_boxplot(width=0.1) + 
  xlab("Method") + ylab("Type I Error") + ggtitle("Cutoff at 0.2") + 
  geom_hline(yintercept=0.2, linetype="dashed", color="red", linewidth=1)+ylim(0, 0.3)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=14))
