suppressMessages(library(dirmult))
suppressMessages(library(dplyr))
suppressMessages(library(MASS))

count_permutation <- function(real_count_mat, nSample=100, seed=1){
  # rows are taxa, columns are samples
  if(nSample > ncol(real_count_mat)){
    nSample <- ncol(real_count_mat)
  }
  set.seed(seed)
  permuted_count <- apply(real_count_mat, 1, 
                          function(taxa_counts) taxa_counts[sample.int(ncol(real_count_mat), nSample)]) |> t()
  return(permuted_count)
}

# gen_dirparam_multiply <- function(baseline, propDA=0.1, balanced=T, enrich=F,mean_logfold = 1,
#                                   sd_logfold=0, seed=1){
#   g1_gamma <- baseline
#   g2_gamma <- baseline
#   nTaxa <- length(g1_gamma)
#   DAtaxa <- sample(seq(1, nTaxa), size=round(nTaxa * propDA))
# 
#   avglogfold <- mean_logfold
#   if (enrich){
#     avglogfold <- abs(mean_logfold)
#   } else{
#     avglogfold <- -abs(mean_logfold)
#   }
#   
#   DA_logfoldchanges <- rnorm(round(nTaxa * propDA), 
#                           mean=avglogfold, sd=sd_logfold)
# 
#   if(balanced){
#     subset_DA <- sample(seq(1, length(DA_logfoldchanges)), length(DA_logfoldchanges)/2)
#     DA_logfoldchanges[subset_DA] <- -DA_logfoldchanges[subset_DA]
#   }
#   
#   logfoldchanges <- rep(0, nTaxa)
#   logfoldchanges[DAtaxa] <- DA_logfoldchanges
#   g2_gamma <- g2_gamma * exp(logfoldchanges)
#   
#   
#   result <- list(g1_gamma=g1_gamma, g2_gamma=g2_gamma, log_fold_change = logfoldchanges)
#   return(result)
#   
# }


SimulateCount <- function (
  dirmult_composition, nSample = 100, seqdepth_mean=20000, seqdepth_theta=10, depth_confound=F,
  propDA=0.1, mean_logfold = log(4), sd_logfold=0.1, direction =c("balanced", "unbalanced"),
  zero_inflation=0, covar_confound=F, seed=1
){
  # check some parameters
  DA_direction <- match.arg(direction)
  stopifnot(zero_inflation < 1 | zero_inflation >= 0)
  stopifnot(propDA <= 0.5 | propDA >= 0)
  # simulate count based on given baseline composition parameter and sequencing depth parameters
  nTaxa <- length(dirmult_composition)
  
  
  # separate the samples into two groups
  grp_proportion <- 0.5
  g1_samplesize <- round(nSample * grp_proportion)
  g2_samplesize <- nSample - g1_samplesize
  
  # create a matrix in which each row is the dirichlet distribution parameter for a sample composition
  all_dirparams <- replicate(nSample, dirmult_composition) |> t()

  ## we multiply the Dirichlet parameters with fold changes determined by exposure variable and confounding variable
  logfolds <- matrix(0, nrow=nSample, ncol=nTaxa)
  
  ## main exposure
  main_exposure_vec <- c(rep(0, g1_samplesize), rep(1, g2_samplesize)) # binary covariate of each sample
  main_exposure_eff_vec <- rep(0, nTaxa) # effect sizes of each taxon
  if (propDA > 0){
    main_exposure_mat <- replicate(nTaxa, main_exposure_vec)
    n_DAtaxa <- round(nTaxa * propDA)
    selected_DAtaxa <- sample(nTaxa, n_DAtaxa)
    DA_effsizes <- rnorm(n_DAtaxa, mean=mean_logfold, sd=sd_logfold)
    if (DA_direction == "balanced"){
      mask <- rep(1, n_DAtaxa)
      mask[sample(n_DAtaxa, n_DAtaxa/2)] <- -1
      DA_effsizes <- DA_effsizes * mask
    }
    main_exposure_eff_vec[selected_DAtaxa] <- DA_effsizes
    main_exposure_eff_mat <- replicate(nSample, main_exposure_eff_vec) |> t()
    logfolds <- logfolds + main_exposure_mat * main_exposure_eff_mat
  }
  
  
  ## TODO: add confounding variables
  
  
  ## sum up the contributions from main exposure and confounding variable
  all_dirparams <- exp(logfolds) * all_dirparams

  
  sim_metadata <- data.frame(exposure = main_exposure_vec) # TODO: confounding factor
  rownames(sim_metadata) <- sprintf("Indv%d", seq(1, nSample))
  sim_taxainfo <- data.frame(exposure_eff = main_exposure_eff_vec) # TODO: confounding factor
  rownames(sim_taxainfo) <- sprintf("ASV%d", seq(1, nTaxa))
  
  set.seed(seed)
  all_seqdepths <- rep(0, nSample)
  # add seqdepth confounding if required
  if (!depth_confound){
    all_seqdepths <- rnegbin(n=nSample, mu=seqdepth_mean, theta=seqdepth_theta) 
  } else{
    seqdepths_g1 <- rnegbin(n=g1_samplesize, mu=seqdepth_mean/2, theta=seqdepth_theta)
    seqdepths_g2 <- rnegbin(n=g2_samplesize, mu=seqdepth_mean*2, theta=seqdepth_theta)
    all_seqdepths <- c(seqdepths_g1, seqdepths_g2)
  }
  # simulate count matrix
  sim_otu_table <- matrix(0, nrow=nSample, ncol=nTaxa)
  for (j in 1:nSample){
    composition <- rdirichlet(n=1, alpha=all_dirparams[j, ]) |> as.vector()
    sim_otu_table[j, ] <- rmultinom(1, size=all_seqdepths[j], prob=composition)
  }

  # add zero inflation if required
  if (zero_inflation > 0){
    nonzero_indices <- which(sim_otu_table > 0)
    truncate_indices <- sample(nonzero_indices, length(nonzero_indices) * zero_inflation)
    sim_otu_table[truncate_indices] <- 0
  }
  
  colnames(sim_otu_table) <- sprintf("ASV%d", seq(1, nTaxa))
  rownames(sim_otu_table) <- sprintf("Indv%d", seq(1, nSample))
  
  output <- list(count_mat = sim_otu_table,
                 sample_metadata = sim_metadata,
                 taxa_info = sim_taxainfo)
  return(output)
}



# SimCount2 <- function (
#     dirmult_composition, seqdepth_mean, seqdepth_theta=10, nSample = 100,
#     grp_ratio = 1, propDA=0.1, mean_logfold = 1,
#     sd_logfold=0, seed=1,
#     DA_direction=c("enrich", "deplete", "mix")
# ){
#   
#   set.seed(seed)
#   DA_direction <- match.arg(DA_direction)
#   if(propDA > 0.5){
#     stop("The number of DA taxa should not exceed 50% of all the taxa!")
#   }
#   # simulate count based on given baseline composition parameter and sequencing depth parameters
#   nTaxa <- length(dirmult_composition)
#   grp_proportion <- grp_ratio / (grp_ratio + 1) # the sample proportion in two constrasting groups
#   
#   # separate the samples into two groups
#   g1_samplesize <- round(nSample * grp_proportion)
#   g2_samplesize <- nSample - g1_samplesize
#   ## covariate vector
#   covariate_vec <- c(rep(0, g1_samplesize), rep(1, g2_samplesize))
#   ## TODO: confounder to be added?
#   
#   # generate compositions for all the samples
#   DAtaxa <- sample(seq(1, nTaxa), size=round(nTaxa * propDA))
#   ## generate fold changes
#   avglogfold <- mean_logfold
#   if (DA_direction == "enrich"){
#     avglogfold <- abs(mean_logfold)
#   } else if (DA_direction == "deplete"){
#     avglogfold <- -abs(mean_logfold)
#   }
#   DA_logfoldchanges <- rnorm(round(nTaxa * propDA), 
#                              mean=avglogfold, sd=sd_logfold)
#   if(DA_direction == "mix"){
#     subset_DA <- sample(seq(1, length(DA_logfoldchanges)), length(DA_logfoldchanges)/2)
#     DA_logfoldchanges[subset_DA] <- -DA_logfoldchanges[subset_DA]
#   }
#   logfoldchanges <- rep(0, nTaxa)
#   logfoldchanges[DAtaxa] <- DA_logfoldchanges
#   expanded_foldchanges <- exp(outer(covariate_vec, logfoldchanges))
#   
#   ## generate compositions
#   compositions <- rdirichlet(n=nSample, dirmult_composition)
#   compositions <- compositions * expanded_foldchanges
#   compositions <- compositions / rowSums(compositions)
# 
#   all_seqdepths <- rnegbin(n=nSample, mu=seqdepth_mean, theta=seqdepth_theta)
#   # simulate count matrix
#   sim_otu_table <- matrix(0, nrow=nTaxa, ncol=nSample)
#   for (j in 1:nSample){
#     sim_otu_table[, j] <- rmultinom(1, size=all_seqdepths[j], prob=compositions[j, ])
#   }
#   
#   sim_metadata <- data.frame(Group=covariate_vec)
#   rownames(sim_metadata) <- sprintf("Indv%d", seq(1, nSample))
#   sim_taxainfo <- data.frame(Logfold=logfoldchanges)
#   rownames(sim_taxainfo) <- sprintf("OTU%d", seq(1, nTaxa))
#   
#   rownames(sim_otu_table) <- sprintf("OTU%d", seq(1, nTaxa))
#   colnames(sim_otu_table) <- sprintf("Indv%d", seq(1, nSample))
#   
#   output <- list(count_mat = sim_otu_table,
#                  sample_metadata = sim_metadata,
#                  taxa_info = sim_taxainfo)
#   return(output)
# }

