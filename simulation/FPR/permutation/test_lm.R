rm(list=ls())
library(POLDA)
# library(dacomp)
library(survival)

real_data_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data"
simulation_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation"
source(file.path(simulation_folder, "FPR", "simulation_utils.R"))
global_gut_infant <- readRDS(file.path(real_data_folder, "global_gut", "global_gut_infant.rds"))
global_gut_adult <- readRDS(file.path(real_data_folder, "global_gut", "global_gut_adult.rds"))


count_matrix <- global_gut_adult$otu_table
metadata <- global_gut_adult$metadata
metadata$isUSA <- metadata$country == "USA"

count_USA <- count_matrix[, metadata$isUSA]
prevalences_USA <- rowMeans(count_USA != 0)
count_USA <- count_USA[prevalences_USA > 0.2, ]


metadata <- data.frame(X=rbinom(129, 1, 0.5))
count_permutation <- function(real_count_mat, seed=1){
  # rows are taxa, columns are samples
  
  nSample <- ncol(real_count_mat)
  set.seed(seed)
  permuted_count <- apply(real_count_mat, 1,
                          function(taxa_counts) taxa_counts[sample.int(ncol(real_count_mat), nSample)]) |> t()
  return(permuted_count)
}

j <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
permuted_counts <- count_permutation(count_USA, seed=j)
# permuted_counts <- truncate_counts(permuted_counts, targetzeros = 0.7)
prevalences <- rowMeans(permuted_counts > 0)
permuted_counts <- permuted_counts[prevalences > 0.1, ]

# selected_DA_taxa <- polda_pvals$Taxon[which.min(polda_pvals$pval)]
refcounts <- colSums(permuted_counts)


lm_pvals <- rep(0, nrow(permuted_counts))


for (k in 1:nrow(permuted_counts)){

  count_indv <- permuted_counts[k,]
  logratios <- rep(0, length(refcounts))
  logratios[count_indv != 0] <- log(count_indv[count_indv!=0] / refcounts[count_indv!=0])
  logratios[count_indv == 0] <- log(1 / (refcounts[count_indv==0] + 1))

  covariate <- metadata$X

  linear_model <- lm(logratios ~ metadata$X) |> summary()
  lm_pvals[k] <- linear_model$coefficients[2,4]
  
}


avg_5 <- mean(lm_pvals < 0.05)
avg_10 <- mean(lm_pvals < 0.1)
avg_20 <- mean(lm_pvals < 0.2)
avg_30 <- mean(lm_pvals < 0.3)

output <- sprintf("%d\t%.4f\t%.4f\t%.4f\t%.4f", j, avg_5, avg_10, avg_20, avg_30)
write(output, file=file.path(simulation_folder, "FPR", "permutation", "lm_FPR.txt"),append=T)


