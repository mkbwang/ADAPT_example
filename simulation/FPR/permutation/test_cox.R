rm(list=ls())
library(POLDA)
# library(dacomp)
library(survival)

real_data_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data"
simulation_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation"
source(file.path(simulation_folder, "FPR", "simulation_utils.R"))
global_gut_adult <- readRDS(file.path(real_data_folder, "global_gut", "global_gut_USadult.rds"))


count_matrix <- global_gut_adult$otu_table
metadata <- global_gut_adult$metadata
# metadata$isUSA <- metadata$country == "USA"
# 
# count_USA <- count_matrix[, metadata$isUSA]
# prevalences_USA <- rowMeans(count_USA != 0)
# count_USA <- count_USA[prevalences_USA > 0.05, ] |> as.matrix()
# count_USA <- as.matrix(count_USA)


j <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(j)
metadata$X <- sample(metadata$X, nrow(metadata))
# count_permutation <- function(real_count_mat, seed=1){
#   # rows are taxa, columns are samples
# 
#   nSample <- ncol(real_count_mat)
#   set.seed(seed)
#   permuted_count <- apply(real_count_mat, 1,
#                           function(taxa_counts) taxa_counts[sample.int(ncol(real_count_mat), nSample)]) |> t()
#   return(permuted_count)
# }
# 
# j <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# permuted_counts <- count_permutation(count_USA, seed=j)
# # permuted_counts <- truncate_counts(permuted_counts, targetzeros = 0.7)
# prevalences <- rowMeans(permuted_counts > 0)
# permuted_counts <- permuted_counts[prevalences > 0.1, ]

refcounts <- colSums(count_matrix)

cox_pvals <- rep(0, nrow(count_matrix))

for (k in 1:nrow(count_matrix)){

  count_indv <- count_matrix[k,]
  logratios <- rep(0, length(refcounts))
  logratios[count_indv != 0] <- log(count_indv[count_indv!=0] / refcounts[count_indv!=0])
  logratios[count_indv == 0] <- log(1 / (refcounts[count_indv==0] + 1))

  covariate <- metadata$X

  cox_model <- coxph(Surv(exp(-logratios), count_indv > 0) ~ covariate) |> summary()
  cox_pvals[k] <- cox_model$coefficients[1,5]
  
}


avg_5 <- mean(cox_pvals < 0.05)
avg_10 <- mean(cox_pvals < 0.1)
avg_20 <- mean(cox_pvals < 0.2)
avg_30 <- mean(cox_pvals < 0.3)
output <- sprintf("%d\t%.4f\t%.4f\t%.4f\t%.4f", j, avg_5, avg_10, avg_20, avg_30)


write(output, file=file.path(simulation_folder, "FPR", "permutation", "cox_FPR.txt"),append=T)

