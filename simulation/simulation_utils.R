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


SimulateCount <- function (
  baseline_dir, nSample = 100, nTaxa=length(baseline_dir), seqdepth_mean=20000, seqdepth_theta=10, depth_confound=F,
  propDA=0.1, main_mean_logfold = log(4), main_sd_logfold=0.1, direction =c("balanced", "unbalanced"),
  propcf = 0, cf_sd=1, cf_main_corr=0.5, cf_mean_logfold=1, cf_sd_logfold=0.1,
  zero_inflation=0, seed=1
){
  # check parameters validity
  DA_direction <- match.arg(direction)
  stopifnot(zero_inflation < 1 & zero_inflation >= 0)
  stopifnot(propDA < 0.5 & propDA >= 0)
  stopifnot(cf_main_corr >= 0 & cf_main_corr < 1)
  stopifnot(nTaxa <= length(baseline_dir))

  template_dim <- length(baseline_dir)
  # create a matrix in which each row is the dirichlet distribution parameter for a sample composition
  all_dirparams <- replicate(nSample, baseline_dir) |> t()

  ## we multiply the Dirichlet parameters with fold changes determined by exposure variable and confounding variable
  logfolds <- matrix(0, nrow=nSample, ncol=template_dim)
  
  ## first generate the continuous and potential confounder
  cf_vec <- rnorm(nSample, mean=0, sd=cf_sd)
  ## then generate the binary main variable using an intermediate variable
  intermediate <- sqrt(cf_main_corr) * cf_vec + sqrt(1 - cf_main_corr) * rnorm(nSample, mean=0, sd=cf_sd)
  main_vec <- rep(0, nSample)
  main_vec[intermediate > median(intermediate)] <- 1

  
  main_eff_vec <- rep(0, template_dim) # main effect sizes of each taxon
  if (propDA > 0){ # exists at least one taxa whose abundance is related to the main binary covariate
    n_DAtaxa <- round(template_dim * propDA)
    selected_DAtaxa <- sample(template_dim, n_DAtaxa)
    DA_effsizes <- rnorm(n_DAtaxa, mean=main_mean_logfold, sd=main_sd_logfold)
    if (DA_direction == "balanced"){
      mask <- rep(1, n_DAtaxa)
      mask[sample(n_DAtaxa, n_DAtaxa/2)] <- -1
      DA_effsizes <- DA_effsizes * mask
    }
    main_eff_vec[selected_DAtaxa] <- DA_effsizes
    logfolds <- logfolds + outer(main_vec, main_eff_vec)
  }
  
  cf_eff_vec <- rep(0, template_dim) # confounding effect sizes of each taxon
  if (propcf > 0){ # exists at least one taxa whose abundance is related to the confounder
    n_cftaxa <- round(template_dim * propcf)
    selected_cftaxa <- sample(template_dim, n_cftaxa)
    cf_effsizes <- rnorm(n_cftaxa, mean=cf_mean_logfold, sd=cf_sd_logfold)
    negate_indices <- sample(n_cftaxa, size=round(n_cftaxa/2), replace=F)
    cf_effsizes[negate_indices] <- -cf_effsizes[negate_indices]
    cf_eff_vec[selected_cftaxa] <- cf_effsizes
    logfolds <- logfolds + outer(cf_vec, cf_eff_vec)
  }
  
  
  ## sum up the contributions from main exposure and confounding variable
  all_dirparams <- exp(logfolds) * all_dirparams

  sim_metadata <- data.frame(main = main_vec, confounder=cf_vec) 
  rownames(sim_metadata) <- sprintf("Indv%d", seq(1, nSample))
  sim_taxainfo <- data.frame(log_main_eff = main_eff_vec, log_cf_eff = cf_eff_vec)
  sim_taxainfo$isDA <- main_eff_vec != 0
  sim_taxainfo$iscf <- cf_eff_vec != 0
  sim_taxainfo$taxaname <- sprintf("Taxon%d", seq(1, template_dim))
  rownames(sim_taxainfo) <- sprintf("Taxon%d", seq(1, template_dim))
  
  set.seed(seed)
  all_seqdepths <- rep(0, nSample)
  # add seqdepth confounding if required
  if (!depth_confound){
    all_seqdepths <- rnegbin(n=nSample, mu=seqdepth_mean, theta=seqdepth_theta) 
  } else{
    all_seqdepths[main_vec == 0] <- rnegbin(n=sum(main_vec == 0), mu=seqdepth_mean/2, theta=seqdepth_theta)
    all_seqdepths[main_vec == 1] <- rnegbin(n=sum(main_vec == 1), mu=seqdepth_mean*2, theta=seqdepth_theta)
  }
  # simulate count matrix
  sim_otu_table <- matrix(0, nrow=nSample, ncol=template_dim)
  for (j in 1:nSample){
    composition <- rdirichlet(n=1, alpha=all_dirparams[j, ]) |> as.vector()
    counts_vec <- rmultinom(1, size=all_seqdepths[j], prob=composition) |> as.vector()

    if (zero_inflation > 0){ # truncate some compositions to be exactly zero
      nonzero_entries <- which(counts_vec > 0)
      num_nonzero <- round(length(nonzero_entries) * (1-zero_inflation))
      subset_nonzero_entries <- sample(nonzero_entries, 
                                num_nonzero, 
                                replace=F, 
                                prob=composition[nonzero_entries])
      zero_entries <- setdiff(seq(1, template_dim), subset_nonzero_entries)
      counts_vec[zero_entries] <- 0
    }
    sim_otu_table[j, ] <- counts_vec
  }
  
  colnames(sim_otu_table) <- sprintf("Taxon%d", seq(1, template_dim))
  rownames(sim_otu_table) <- sprintf("Indv%d", seq(1, nSample))
  
  if (nTaxa < length(baseline_dir)){ # take a subset of taxa
    sim_taxainfo <- sim_taxainfo %>% group_by(isDA, iscf) %>% 
      slice_sample(prop=nTaxa/template_dim) %>% as.data.frame()
    rownames(sim_taxainfo) <- sim_taxainfo$taxaname
    sim_otu_table <- sim_otu_table[, sim_taxainfo$taxaname]
  }
  
  output <- list(count_mat = sim_otu_table,
                 sample_metadata = sim_metadata,
                 taxa_info = sim_taxainfo)
  return(output)
}


