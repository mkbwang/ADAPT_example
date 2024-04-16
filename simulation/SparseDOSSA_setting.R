library(SparseDOSSA2)
# get the pretrained template
# original_simulation <- SparseDOSSA2(template = "Stool",  # choose from "Stool", "Vaginal" or "IBD"
#                                  new_features = TRUE,  # should new features be simulated
#                                  n_sample = 100,  # number of samples to simulate
#                                  n_feature = 1000,  # number of features to simulate (when 'new_features = TRUE')
#                                  spike_metadata = "none",
#                                  verbose = F)  # return detailed info
# 
# # filter out abundant taxa (prevalence bigger than 50% and median relative abundance bigger than 2e-4)
# original_template <- original_simulation$template
# original_Ffit <- original_template$F_fit
# original_depth <- original_template$depth_fit
# original_EMfit <- original_template$EM_fit
# original_pi0 <- original_EMfit$fit$pi0 |> unname()
# original_mu <- original_EMfit$fit$mu |> unname()
# original_sigma <- original_EMfit$fit$sigma |> unname()
# 
# 
# original_median_abundance <- rep(0, length(original_pi0))
# nonzero_mask <- original_pi0 < 0.5
# truncated_quantile <- qnorm(0.5 - original_pi0[nonzero_mask])
# original_median_abundance[nonzero_mask] <- exp(original_mu[nonzero_mask] +
#                                                  truncated_quantile*original_sigma[nonzero_mask])
# original_median_relabd <- original_median_abundance / sum(original_median_abundance)
# subset_taxa <- which(original_median_relabd > 2e-4)
# 
# 
# subset_pi0 <- original_pi0[subset_taxa]
# subset_mu <- original_mu[subset_taxa]
# subset_sigma <- original_sigma[subset_taxa]
# 
# subset_EMfit <- list(fit=list(pi0=subset_pi0, 
#                               mu=subset_mu, 
#                               sigma=subset_sigma))
# 
# subset_template <- list(EM_fit=subset_EMfit,
#                         F_fit=original_Ffit,
#                         depth_fit=original_depth)
# saveRDS(subset_template, file=file.path("simulation/sparsedossa_baseline.rds"))
# 
# Stool_simulation <- SparseDOSSA2(template = subset_template,  # choose from "Stool", "Vaginal" or "IBD"
#                                  new_features = TRUE,  # should new features be simulated
#                                  n_sample = 50,  # number of samples to simulate
#                                  n_feature = 1000,  # number of features to simulate (when 'new_features = TRUE')
#                                  spike_metadata = "none",
#                                  median_read_depth=5e4,
#                                  verbose = F)  # return detailed info
# 
# 
# simulated_countmat <- Stool_simulation$simulated_data
# prevalences <- rowMeans(simulated_countmat > 0)


baseline_params <- readRDS('/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation/sparsedossa_baseline.rds')


# simulation settings
SimulateSettings <- function (baseline,
    nSample = 100, nTaxa=1000, seqdepth_logmean=log(5e4), seqdepth_logsd=1,
    propDA=0.1, main_mean_logfold = log(4), main_sd_logfold=0.1, direction =c("balanced", "unbalanced"),
    propcf = 0, cf_sd=1, cf_main_corr=0.5, cf_mean_logfold=1, cf_sd_logfold=0.1, seed=1
){
  
  # check parameters validity
  DA_direction <- match.arg(direction)
  stopifnot(propDA < 0.9 & propDA >= 0)
  stopifnot(cf_main_corr >= 0 & cf_main_corr < 1)

  set.seed(seed)
  # spike_main_metadata <- data.frame(metadata_datum=integer(), 
  #                        feature_spiked=character(), 
  #                        associated_property=character(), 
  #                        effect_size=numeric())
  
  ## first generate the continuous and potential confounder
  cf_vec <- rnorm(nSample, mean=0, sd=cf_sd)
  ## then generate the binary main variable using an intermediate variable
  intermediate <- sqrt(cf_main_corr) * cf_vec + sqrt(1 - cf_main_corr) * rnorm(nSample, mean=0, sd=cf_sd)
  main_vec <- rep(0, nSample)
  main_vec[intermediate > median(intermediate)] <- 1
  
  metadata_mat <- cbind(main_vec, cf_vec)
  
  ## effect sizes of the main covariate
  spike_main_metadata <- NULL
  
  main_eff_vec <- rep(0, nTaxa) # main effect sizes of each taxon
  if (propDA > 0){ # exists at least one taxa whose abundance is related to the main binary covariate
    n_DAtaxa <- round(nTaxa * propDA)
    selected_DAtaxa <- sample(nTaxa, n_DAtaxa)
    DA_effsizes <- rnorm(n_DAtaxa, mean=main_mean_logfold, sd=main_sd_logfold)
    if (DA_direction == "balanced"){
      mask <- rep(1, n_DAtaxa)
      mask[sample(n_DAtaxa, n_DAtaxa/2)] <- -1
      DA_effsizes <- DA_effsizes * mask
    }
    main_eff_vec[selected_DAtaxa] <- DA_effsizes
    spike_main_metadata <- data.frame(metadata_datum = rep(1, n_DAtaxa*2),
                                      feature_spiked = sprintf("Feature%d", rep(selected_DAtaxa, 2)),
                                      associated_property = rep(c("abundance", "prevalence"), each=n_DAtaxa),
                                      effect_size = rep(DA_effsizes, 2))
  }
  
  ## effect sizes of the confounding covariate
  spike_cf_metadata <- NULL
  cf_eff_vec <- rep(0, nTaxa) # confounding effect sizes of each taxon
  if (propcf > 0){ # exists at least one taxa whose abundance is related to the confounder
    n_cftaxa <- round(nTaxa * propcf)
    selected_cftaxa <- sample(nTaxa, n_cftaxa)
    cf_effsizes <- rnorm(n_cftaxa, mean=cf_mean_logfold, sd=cf_sd_logfold)
    negate_indices <- sample(n_cftaxa, size=round(n_cftaxa/2), replace=F)
    cf_effsizes[negate_indices] <- -cf_effsizes[negate_indices]
    cf_eff_vec[selected_cftaxa] <- cf_effsizes
    spike_cf_metadata <- data.frame(metadata_datum = rep(2, n_cftaxa*2),
                                      feature_spiked = sprintf("Feature%d", rep(selected_cftaxa, 2)),
                                      associated_property = rep(c("abundance", "prevalence"), each=n_cftaxa),
                                      effect_size = rep(cf_effsizes, 2))
  }
  
  spike_metadata_df <- data.frame(metadata_datum=integer(), 
                                  feature_spiked=character(), 
                                  associated_property=character(), 
                                  effect_size=numeric())
  
  
  if (!is.null(spike_main_metadata)){
    spike_metadata_df <- rbind(spike_metadata_df, spike_main_metadata)
  }
 
  if (!is.null(spike_cf_metadata)){
    spike_metadata_df <- rbind(spike_metadata_df, spike_cf_metadata)
  }
  
  if (nrow(spike_metadata_df) == 0){
    spike_metadata_df <- "none"
  }
  
  # set up sequencing depths parameters
  baseline$depth_fit[1] <- seqdepth_logmean
  baseline$depth_fit[2] <- seqdepth_logsd
  
  
  sim_taxainfo <- data.frame(log_main_eff = main_eff_vec, log_cf_eff = cf_eff_vec)
  sim_taxainfo$isDA <- main_eff_vec != 0
  sim_taxainfo$iscf <- cf_eff_vec != 0
  sim_taxainfo$taxaname <- sprintf("Taxon%d", seq(1, nTaxa))
  rownames(sim_taxainfo) <- sprintf("Taxon%d", seq(1, nTaxa))
  
  
  output <- list(sample_metadata = metadata_mat,
                 taxa_info = sim_taxainfo,
                 baseline_param=baseline,
                 spike_metadata=spike_metadata_df)
  return(output)
}

