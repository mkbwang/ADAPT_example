# example histogram of p values

rm(list=ls())
library(ggplot2)


data_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA/simulation/data'
polda_folder <- '/home/wangmk/UM/Research/MDAWG/POLDA/simulation/POLDA'

source("/home/wangmk/UM/Research/MDAWG/POLDA/POLDA/overdisperse_GLM.R")

good_null_example <- c()
for (j in 1:100){
  null_example <- readRDS(file.path(polda_folder, 
                                    "null", 
                                    sprintf("summary_null_%d.rds", j)))
  if (length(null_example$polda_result$Reference_Taxa) == 1000){
    good_null_example <- c(good_null_example, j)
  }
}

null_example <- readRDS(file.path(polda_folder, 
                                  "null", 
                                  sprintf("summary_null_%d.rds", 21)))
pval_null <- null_example$polda_result$P_Value
dist_null <- data.frame(x=seq(0,1,0.01), y=rep(1, 101))
# hist(pval_null$pval, nclass=30)
pval_null_hist <- ggplot(pval_null, aes(x=pval, y=..density..)) + 
  geom_histogram(breaks=seq(0, 1, 0.05), color="black", fill="white")+
  xlab("Uniform(0,1) Distribution") +ylab("P Value Density") +
  geom_line(data=dist_null, aes(x=x, y=y), linetype="dashed", color = "blue1", size=1.2)+
  ggtitle("No Differentially Abundant Taxa")+theme_bw()+ ylim(0, 2.5)+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5))


DA_example <- readRDS(file.path(data_folder, 
                                "DA_abundant", 
                                "low",
                                "AGP_simulation_unbalanced_abundant_low_20.rds"))

count_DA <- DA_example$count_mat
metadata_DA <- DA_example$sample_metadata
pval_DA <- reference_GLM(count_data=count_DA,
                           metadata=metadata_DA,
                           covar="X")
bumfit <- Bum(pval_DA$pval)
lambda_hat <- bumfit@lhat
a_hat <- bumfit@ahat
x_vec <- seq(0.01, 1, 0.01)
uniform_component <- rep(lambda_hat, length(x_vec))
dist_density <- uniform_component + (1-lambda_hat) * a_hat * x_vec^(a_hat - 1)
dist_df <- data.frame(x = x_vec,
                      y_uniform=uniform_component,
                      y_dist=dist_density)

pval_DA_hist <- ggplot(pval_DA, aes(x=pval, y=..density..)) + 
  geom_histogram(breaks=seq(0, 1, 0.05), color="black", fill="white")+
  xlab("BUM(a=0.192, \u03BB = 0.846) Distribution") +ylab("P Value Density") +
  geom_line(data=dist_df, aes(x=x_vec, y=y_uniform), linetype="dashed", color = "blue1", size=1.2)+
  geom_line(data=dist_df, aes(x=x_vec, y=y_dist), linetype="dashed", color = "chocolate", size=1.2)+
  ggtitle("5% Differentially Abundant Taxa")+theme_bw()+
  theme(text=element_text(size=11),
        plot.title = element_text(hjust = 0.5))

library(cowplot)
combined_plot <- plot_grid(pval_null_hist, pval_DA_hist, ncol=2)

ggsave(filename = "simulation/plots/pval_hist_combined.pdf", 
       plot=combined_plot,
       device=cairo_pdf, width=9.06, height=3.66, units="in",
       dpi="print")
