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


## generate dirichlet distribution for two groups
# gen_dirparam_swap <- function(baseline, propDA=0.1, balanced=T, enrich=F, seed=1){
#   stopifnot(propDA <= 0.2)
#   num_taxa <- length(baseline)
#   num_DAtaxa <- round(num_taxa * propDA)
#   param_20qt <- quantile(baseline, probs=0.2)
#   param_80qt <- quantile(baseline, probs=0.8)
#   rare_taxa_ID <- which(baseline < param_20qt)
#   prevalent_taxa_ID <- which(baseline > param_80qt)
#   
#   g1_gamma <- baseline
#   g2_gamma <- baseline
#   set.seed(seed)
#   if(balanced){ # some rare taxa enriched and some prevalent taxa depleted
#     DA_raretaxa_ID <- sample(rare_taxa_ID, size=num_DAtaxa/2)
#     raretaxa_gamma <- gamma_param[DA_raretaxa_ID]
#     DA_prevalenttaxa_ID <- sample(prevalent_taxa_ID, size=num_DAtaxa/2)
#     prevalenttaxa_gamma <- gamma_param[DA_prevalenttaxa_ID]
#     g2_gamma[DA_raretaxa_ID] <- sample(prevalenttaxa_gamma, size=num_DAtaxa/2)
#     g2_gamma[DA_prevalenttaxa_ID] <- sample(raretaxa_gamma, size=num_DAtaxa/2)
#   } else{# either rare taxa enriched or prevalent taxa depleted, one direction
#     DA_raretaxa_ID <- sample(rare_taxa_ID, size=num_DAtaxa)
#     raretaxa_gamma <- gamma_param[DA_raretaxa_ID]
#     DA_prevalenttaxa_ID <- sample(prevalent_taxa_ID, size=num_DAtaxa)
#     prevalenttaxa_gamma <- gamma_param[DA_prevalenttaxa_ID]
#     if (enrich){ # rare taxa get enriched
#       g2_gamma[DA_raretaxa_ID] <- sample(prevalenttaxa_gamma, size=num_DAtaxa)
#     } else{ # prevalent taxa get depleted
#       g2_gamma[DA_prevalenttaxa_ID] <- sample(raretaxa_gamma, size=num_DAtaxa)
#     }
#   }
#   result <- list(g1_gamma=g1_gamma, g2_gamma=g2_gamma, log_fold_change = log(g2_gamma/g1_gamma))
#   return(result)
# }

gen_dirparam_multiply <- function(baseline, propDA=0.1, balanced=T, enrich=F,mean_logfold = 1,
                                  sd_logfold=0, seed=1){
  g1_gamma <- baseline
  g2_gamma <- baseline
  nTaxa <- length(g1_gamma)
  DAtaxa <- sample(seq(1, nTaxa), size=round(nTaxa * propDA))

  avglogfold <- mean_logfold
  if (enrich){
    avglogfold <- abs(mean_logfold)
  } else{
    avglogfold <- -abs(mean_logfold)
  }
  
  DA_logfoldchanges <- rnorm(round(nTaxa * propDA), 
                          mean=avglogfold, sd=sd_logfold)

  if(balanced){
    subset_DA <- sample(seq(1, length(DA_logfoldchanges)), length(DA_logfoldchanges)/2)
    DA_logfoldchanges[subset_DA] <- -DA_logfoldchanges[subset_DA]
  }
  
  logfoldchanges <- rep(0, nTaxa)
  logfoldchanges[DAtaxa] <- DA_logfoldchanges
  g2_gamma <- g2_gamma * exp(logfoldchanges)
  
  
  result <- list(g1_gamma=g1_gamma, g2_gamma=g2_gamma, log_fold_change = logfoldchanges)
  return(result)
  
}


SimulateCount <- function (
  dirmult_composition, seqdepth_mean, seqdepth_theta=10, nSample = 100,
  DA_proportion=0.1, mean_logfold = 1, # confounding= 
  sd_logfold=0, seed=1,
  DA_direction=c("balanced", "unbalanced")
){
  DA_direction <- match.arg(DA_direction)
  if(DA_proportion > 0.2){
    stop("The number of DA taxa should not exceed 20% of all the taxa!")
  }
  # simulate count based on given baseline composition parameter and sequencing depth parameters
  nTaxa <- length(dirmult_composition)
  grp_proportion <- 0.5
  
  # separate the samples into two groups
  g1_samplesize <- round(nSample * grp_proportion)
  g2_samplesize <- nSample - g1_samplesize
  
  # create a matrix in which each row is the dirichlet distribution parameter for a sample composition
  all_dirparams <- replicate(nSample, dirmult_composition) |> t()

  ## we multiply the Dirichlet parameters with fold changes determined by exposure variable and confounding variable
  logfolds <- matrix(0, nrow=nSample, ncol=nTaxa)
  
  
  ## main exposure
  main_exposure_vec <- c(rep(0, g1_samplesize), rep(1, g2_samplesize)) # binary covariate of each sample
  main_exposure_mat <- replicate(nTaxa, main_exposure_vec)
  main_exposure_eff_vec <- rep(0, nTaxa) # effect sizes of each taxon
  n_DAtaxa <- round(nTaxa * DA_proportion)
  selected_DAtaxa <- sample(nTaxa, n_DAtaxa)
  DA_effsizes <- rnorm(n_DAtaxa, mean=mean_logfold, sd=sd_logfold)
  if (DA_direction == "balanced"){
    mask <- rep(1, n_DAtaxa)
    mask[sample(n_DAtaxa, n_DAtaxa/2)] <- -1
    DA_effsizes <- DA_effsizes * mask
  }
  main_exposure_eff_vec[selected_DAtaxa] <- DA_effsizes
  main_exposure_eff_mat <- replicate(nSample, main_exposure_eff_vec) |> t()
  
  
  ## TODO: add confounding variables
  
  
  ## sum up the contributions from main exposure and confounding variable
  logfolds <- logfolds + main_exposure_mat * main_exposure_eff_mat
  all_dirparams <- exp(logfolds) * all_dirparams
  
  # sim_dirparams <- NULL
  # if (DA_direction == "balanced"){
  #   # Are all the DA taxa changing in one direction?
  #   sim_dirparams <- gen_dirparam_multiply(baseline=dirmult_composition, propDA=DA_proportion, mean_logfold = mean_logfold, balanced=T, seed=seed)
  # } else {
  #   sim_dirparams <- gen_dirparam_multiply(baseline=dirmult_composition, propDA=DA_proportion, mean_logfold = mean_logfold, balanced=F, enrich = T, seed=seed)
  # }
  
  # draw compositions from Dirichlet distributions
  # set.seed(seed)
  # g1_compositions <- rdirichlet(n=g1_samplesize, sim_dirparams$g1_gamma) |> t()
  # g2_compositions <- rdirichlet(n=g2_samplesize, sim_dirparams$g2_gamma) |> t()
  # all_compositions <- cbind(g1_compositions, g2_compositions)
  
  sim_metadata <- data.frame(exposure = main_exposure_vec) # TODO: confounding factor
  rownames(sim_metadata) <- sprintf("Indv%d", seq(1, nSample))
  sim_taxainfo <- data.frame(exposure_eff = main_exposure_eff_vec) # TODO: confounding factor
  rownames(sim_taxainfo) <- sprintf("ASV%d", seq(1, nTaxa))
  
  set.seed(seed)
  ## TODO: seqdepth confounding 
  all_seqdepths <- rnegbin(n=nSample, mu=seqdepth_mean, theta=seqdepth_theta)
  # simulate count matrix
  sim_otu_table <- matrix(0, nrow=nSample, ncol=nTaxa)
  for (j in 1:nSample){
    composition <- rdirichlet(n=1, alpha=all_dirparams[j, ]) |> as.vector()
    sim_otu_table[j, ] <- rmultinom(1, size=all_seqdepths[j], prob=composition)
  }

  colnames(sim_otu_table) <- sprintf("ASV%d", seq(1, nTaxa))
  rownames(sim_otu_table) <- sprintf("Indv%d", seq(1, nSample))
  
  output <- list(count_mat = sim_otu_table,
                 sample_metadata = sim_metadata,
                 taxa_info = sim_taxainfo)
  return(output)
}



SimCount2 <- function (
    dirmult_composition, seqdepth_mean, seqdepth_theta=10, nSample = 100,
    grp_ratio = 1, DA_proportion=0.1, mean_logfold = 1,
    sd_logfold=0, seed=1,
    DA_direction=c("enrich", "deplete", "mix")
){
  
  set.seed(seed)
  DA_direction <- match.arg(DA_direction)
  if(DA_proportion > 0.5){
    stop("The number of DA taxa should not exceed 50% of all the taxa!")
  }
  # simulate count based on given baseline composition parameter and sequencing depth parameters
  nTaxa <- length(dirmult_composition)
  grp_proportion <- grp_ratio / (grp_ratio + 1) # the sample proportion in two constrasting groups
  
  # separate the samples into two groups
  g1_samplesize <- round(nSample * grp_proportion)
  g2_samplesize <- nSample - g1_samplesize
  ## covariate vector
  covariate_vec <- c(rep(0, g1_samplesize), rep(1, g2_samplesize))
  ## TODO: confounder to be added?
  
  # generate compositions for all the samples
  DAtaxa <- sample(seq(1, nTaxa), size=round(nTaxa * DA_proportion))
  ## generate fold changes
  avglogfold <- mean_logfold
  if (DA_direction == "enrich"){
    avglogfold <- abs(mean_logfold)
  } else if (DA_direction == "deplete"){
    avglogfold <- -abs(mean_logfold)
  }
  DA_logfoldchanges <- rnorm(round(nTaxa * DA_proportion), 
                             mean=avglogfold, sd=sd_logfold)
  if(DA_direction == "mix"){
    subset_DA <- sample(seq(1, length(DA_logfoldchanges)), length(DA_logfoldchanges)/2)
    DA_logfoldchanges[subset_DA] <- -DA_logfoldchanges[subset_DA]
  }
  logfoldchanges <- rep(0, nTaxa)
  logfoldchanges[DAtaxa] <- DA_logfoldchanges
  expanded_foldchanges <- exp(outer(covariate_vec, logfoldchanges))
  
  ## generate compositions
  compositions <- rdirichlet(n=nSample, dirmult_composition)
  compositions <- compositions * expanded_foldchanges
  compositions <- compositions / rowSums(compositions)

  all_seqdepths <- rnegbin(n=nSample, mu=seqdepth_mean, theta=seqdepth_theta)
  # simulate count matrix
  sim_otu_table <- matrix(0, nrow=nTaxa, ncol=nSample)
  for (j in 1:nSample){
    sim_otu_table[, j] <- rmultinom(1, size=all_seqdepths[j], prob=compositions[j, ])
  }
  
  sim_metadata <- data.frame(Group=covariate_vec)
  rownames(sim_metadata) <- sprintf("Indv%d", seq(1, nSample))
  sim_taxainfo <- data.frame(Logfold=logfoldchanges)
  rownames(sim_taxainfo) <- sprintf("OTU%d", seq(1, nTaxa))
  
  rownames(sim_otu_table) <- sprintf("OTU%d", seq(1, nTaxa))
  colnames(sim_otu_table) <- sprintf("Indv%d", seq(1, nSample))
  
  output <- list(count_mat = sim_otu_table,
                 sample_metadata = sim_metadata,
                 taxa_info = sim_taxainfo)
  return(output)
}

