library(fitdistrplus)
library(dirmult)
library(MASS)
library(foreach)


# draw from dirichlet distribution parameters
rdirichlet <- function(dir_gamma, nSample, seed=2023){
  set.seed(seed)
  # first generate gamma distribution
  simulation <- sapply(dir_gamma, function(gamma_shape) rgamma(nSample, shape=gamma_shape, rate=1))
  # then normalize to composition
  composition <- apply(simulation, 1, function(sample_count) sample_count/sum(sample_count))
  return(composition)
}

# permute the counts of a real data
count_permutation <- function(real_count_mat, seed=1){
  # rows are taxa, columns are samples
  
  nSample <- ncol(real_count_mat)
  set.seed(seed)
  permuted_count <- apply(real_count_mat, 1,
                          function(taxa_counts) taxa_counts[sample.int(ncol(real_count_mat), nSample)]) |> t()
  return(permuted_count)
}

# estimate dirichlet distribution parameters based on a real data
estimate_param <- function(real_count_mat, seed=1){
  # estimate dirichlet multinomial distribution parameter and library sizes based on real data
  # taxa are rows
  set.seed(seed)
  Nsample <- ncol(real_count_mat)
  Ntaxa <- nrow(real_count_mat)
  permuted_count <- count_permutation(real_count_mat = real_count_mat,
                                      seed=seed)
  
  # negative binomial distribution for sequencing depths
  permuted_depths <- colSums(permuted_count)
  NB_depth_fit <- fitdist(permuted_depths, "nbinom")
  NB_depth_mean <- NB_depth_fit$estimate[2]
  NB_depth_theta <- NB_depth_fit$estimate[1]
  
  # estimate dirichlet multinomial distribution
  dirichlet_param <- dirmult(t(permuted_count))
  # sort the dirichlet parameters
  dirichlet_gamma <- sort(dirichlet_param$gamma, decreasing=TRUE)
  
  result <- list(dirmult_composition=dirichlet_gamma,
                 Seqdepth_mean=NB_depth_mean,
                 Seqdepth_theta=NB_depth_theta)
}


# if a dataset doesn't have enough zeros, truncate some nonzero values to zero
truncate_counts <- function(countmat, targetzeros=0.6){
  
  stopifnot(targetzeros < 0.9)
  cutoff <- quantile(countmat, probs=0.9)
  currentzero_counts <- sum(countmat == 0)
  targetzero_counts <- round(length(countmat) * targetzeros)
  if(targetzero_counts > currentzero_counts){
    nonzero_indices <- which(countmat  < cutoff & countmat > 0)
    truncate_indices <- sample(nonzero_indices, size=targetzero_counts - currentzero_counts)
    countmat[truncate_indices] <- 0
  }
  return(countmat)
  
}

if(!file.exists(file.path(simulation_folder, "FPR", "global_gut_adult_params.rds"))){
  real_data_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data"
  simulation_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example/simulation"
  global_gut_infant <- readRDS(file.path(real_data_folder, "global_gut", "global_gut_infant.rds"))
  global_gut_adult <- readRDS(file.path(real_data_folder, "global_gut", "global_gut_adult.rds"))
  
  
  count_matrix <- global_gut_adult$otu_table
  metadata <- global_gut_adult$metadata
  metadata$isUSA <- metadata$country == "USA"
  
  count_USA <- count_matrix[, metadata$isUSA]
  prevalences_USA <- rowMeans(count_USA != 0)
  count_USA <- count_USA[prevalences_USA > 0.2, ]
  gut_adult_params <- estimate_param(real_count_mat = count_USA, seed=2023)
  saveRDS(gut_adult_params, file.path(simulation_folder, "FPR", "global_gut_adult_params.rds"))
}

