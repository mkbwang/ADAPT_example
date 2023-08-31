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
gen_dirparam <- function(baseline, propDA=0.1, balanced=T, enrich=F, seed=1){
  stopifnot(propDA <= 0.2)
  num_taxa <- length(baseline)
  num_DAtaxa <- round(num_taxa * propDA)
  param_20qt <- quantile(baseline, probs=0.2)
  param_80qt <- quantile(baseline, probs=0.8)
  rare_taxa_ID <- which(baseline < param_20qt)
  prevalent_taxa_ID <- which(baseline > param_80qt)
  
  g1_gamma <- baseline
  g2_gamma <- baseline
  set.seed(seed)
  if(balanced){ # some rare taxa enriched and some prevalent taxa depleted
    DA_raretaxa_ID <- sample(rare_taxa_ID, size=num_DAtaxa/2)
    raretaxa_gamma <- gamma_param[DA_raretaxa_ID]
    DA_prevalenttaxa_ID <- sample(prevalent_taxa_ID, size=num_DAtaxa/2)
    prevalenttaxa_gamma <- gamma_param[DA_prevalenttaxa_ID]
    g2_gamma[DA_raretaxa_ID] <- sample(prevalenttaxa_gamma, size=num_DAtaxa/2)
    g2_gamma[DA_prevalenttaxa_ID] <- sample(raretaxa_gamma, size=num_DAtaxa/2)
  } else{# either rare taxa enriched or prevalent taxa depleted, one direction
    DA_raretaxa_ID <- sample(rare_taxa_ID, size=num_DAtaxa)
    raretaxa_gamma <- gamma_param[DA_raretaxa_ID]
    DA_prevalenttaxa_ID <- sample(prevalent_taxa_ID, size=num_DAtaxa)
    prevalenttaxa_gamma <- gamma_param[DA_prevalenttaxa_ID]
    if (enrich){ # rare taxa get enriched
      g2_gamma[DA_raretaxa_ID] <- sample(prevalenttaxa_gamma, size=num_DAtaxa)
    } else{ # prevalent taxa get depleted
      g2_gamma[DA_prevalenttaxa_ID] <- sample(raretaxa_gamma, size=num_DAtaxa)
    }
  }
  result <- list(g1_gamma=g1_gamma, g2_gamma=g2_gamma, log_fold_change = log(g2_gamma/g1_gamma))
  return(result)
}



SimulateCount <- function (
  dirmult_composition, seqdepth_mean, seqdepth_theta=10, nSample = 100,
  grp_ratio = 1, DA_proportion=0.1, seed=1,
  DA_direction=c("enrich", "deplete", "mix")
){
  DA_direction <- match.arg(DA_direction)
  if(DA_proportion > 0.2){
    stop("The number of DA taxa should not exceed 20% of all the taxa!")
  }
  # simulate count based on given baseline composition parameter and sequencing depth parameters
  nTaxa <- length(dirmult_composition)
  grp_proportion <- grp_ratio / (grp_ratio + 1) # the sample proportion in two constrasting groups
  
  # separate the samples into two groups
  g1_samplesize <- round(nSample * grp_proportion)
  g2_samplesize <- nSample - g1_samplesize
  

  sim_dirparams <- NULL
  if (DA_direction == "mix"){
    # Are all the DA taxa changing in one direction?
    sim_dirparams <- gen_dirparam(baseline=dirmult_composition, propDA=DA_proportion, balanced=T, seed=seed)
  } else if (DA_direction == "enrich"){
    sim_dirparams <- gen_dirparam(baseline=dirmult_composition, propDA=DA_proportion, balanced=F, enrich = T, seed=seed)
  } else{ # deplete
    sim_dirparams <- gen_dirparam(baseline=dirmult_composition, propDA=DA_proportion, balanced=F, enrich = F, seed=seed)
  }
  
  # draw compositions from Dirichlet distributions
  set.seed(seed)
  g1_compositions <- rdirichlet(n=g1_samplesize, sim_dirparams$g1_gamma) |> t()
  g2_compositions <- rdirichlet(n=g2_samplesize, sim_dirparams$g2_gamma) |> t()
  all_compositions <- cbind(g1_compositions, g2_compositions)
  
  
  sim_metadata <- data.frame(Group=c(rep(0, g1_samplesize), rep(1, g2_samplesize)))
  rownames(sim_metadata) <- sprintf("Indv%d", seq(1, nSample))
  sim_taxainfo <- data.frame(Logfold=sim_dirparams$log_fold_change)
  rownames(sim_taxainfo) <- sprintf("OTU%d", seq(1, nTaxa))
  
  
  all_seqdepths <- rnegbin(n=nSample, mu=seqdepth_mean, theta=seqdepth_theta)
  # simulate count matrix
  sim_otu_table <- matrix(0, nrow=nTaxa, ncol=nSample)
  for (j in 1:nSample){
    sim_otu_table[, j] <- rmultinom(1, size=all_seqdepths[j], prob=all_compositions[, j])
  }

  rownames(sim_otu_table) <- sprintf("OTU%d", seq(1, nTaxa))
  colnames(sim_otu_table) <- sprintf("Indv%d", seq(1, nSample))
  
  output <- list(count_mat = sim_otu_table,
                 sample_metadata = sim_metadata,
                 taxa_info = sim_taxainfo)
  return(output)
}

