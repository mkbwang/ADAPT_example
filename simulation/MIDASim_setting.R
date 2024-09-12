
# ## estimate parameters
# data("count.ibd")
# count.ibd.setup_par <- MIDASim.setup(count.ibd, mode = 'parametric')
# 
# # select subset of taxa as simulation template
# prevalence <- colMeans(count.ibd > 0)
# mu_estimate <- count.ibd.setup_par$mu.est
# sigma_estimate <- count.ibd.setup_par$sigma.est
# q_estimate <- count.ibd.setup_par$Q.est
# 
# stats_taxa <- data.frame(prevalence=prevalence,
#                          mu = mu_estimate,
#                          sigma = sigma_estimate,
#                          q = q_estimate,
#                          colnum=seq(1, ncol(count.ibd)))
# 
# library(dplyr)
# ## sigma small (small variance), high mu (high relative abundance), minimum prevalence 0.05
# subset_stats_taxa <- stats_taxa %>% filter(sigma < 0.6) %>% arrange(desc(mu)) %>%
#   filter(prevalence > 0.05)
# ## top 50
# selected_taxa <- subset_stats_taxa$colnum[1:50]
# 
# 
# 
# # prepare template
# selected_mu <- count.ibd.setup_par$mu.est[selected_taxa]
# selected_sigma <- count.ibd.setup_par$sigma.est[selected_taxa]
# selected_Q <- count.ibd.setup_par$Q.est[selected_taxa]
# selected_tetra_corr <- count.ibd.setup_par$tetra.corr[selected_taxa, selected_taxa]
# selected_rel_corr <- count.ibd.setup_par$corr.rel.corrected[selected_taxa, selected_taxa]
# 
# template <- list(mu = selected_mu, sigma = selected_sigma, Q = selected_Q,
#                  tetra_corr = selected_tetra_corr, rel_corr = selected_rel_corr)
# 
# saveRDS(template, "simulation/midasim_baseline.rds")

library(MIDASim)


# I made a customized function to generate simulated counts for a uniform population to replace MIDASim function
draw_uniform_pop <- function(template_mu, template_sigma, template_Q, libsizes, 
                            template_tetra_corr, template_rel_corr){
  
  num_samples <- length(libsizes)
  num_taxa <- length(template_mu)
  
  prob01_mat <- MIDASim:::get.prob01.mat(template_mu,
                                                   template_sigma,
                                                   template_Q,
                                                   libsizes)
  
  modified_setup <- list(mode="parametric",
                         n.taxa = num_taxa,
                         n.sample=num_samples,
                         prob01.mat=prob01_mat,
                         tetra.corr=template_tetra_corr,
                         corr.rel.corrected=template_rel_corr,
                         mu.est=template_mu,
                         sigma.est=template_sigma,
                         Q.est=template_Q,
                         lib.size=libsizes)
  # get relative abundance
  simulated_relabd <- MIDASim(modified_setup, only.rel=TRUE)
  
  # draw from multinomial distribution and get the final counts
  count_mat <- round(libsizes * simulated_relabd$sim_rel)
  # for (j in 1:num_samples){
  #   count_mat[j, ] <- rmultinom(n=1, size=libsizes[j], prob=simulated_data$sim_rel[j, ])
  # }

  return(count_mat)
  
}


# this function only considers differentially abundant taxa wrt a binary condition. No confounders are involved

midasim_simulate <- function(midasim_setting){
  
  metadata_df <- as.data.frame(midasim_setting$sample_metadata)
  num_sample <- nrow(metadata_df)
  taxa_info <- midasim_setting$taxa_info
  num_taxa <- nrow(taxa_info)
  baseline <- midasim_setting$baseline_param
  
  repeats <- num_taxa / length(baseline$mu)
  template_group1_mu <- rep(baseline$mu, repeats)
  template_group2_mu <- template_group1_mu + taxa_info$log_main_eff
  template_sigma <- rep(baseline$sigma, repeats)
  template_Q <- rep(baseline$Q, repeats)
  template_tetra_corr <- diag(num_taxa)
  template_rel_corr <- diag(num_taxa)
  for (j in 1:repeats){
    start_index <- (j-1)*50+1
    end_index <- j*50
    template_tetra_corr[seq(start_index, end_index), seq(start_index, end_index)] <- baseline$tetra_corr
    template_rel_corr[seq(start_index, end_index), seq(start_index, end_index)] <- baseline$rel_corr
  }
  
  num_group1 <- sum(metadata_df$main_vec == 0)
  libsizes_1 <- exp(rnorm(num_group1, mean=baseline$depth_fit[1], sd=baseline$depth_fit[2])) |> round()
  
  num_group2 <- sum(metadata_df$main_vec == 1)
  libsizes_2 <- exp(rnorm(num_group2, mean=baseline$depth_fit[1], sd=baseline$depth_fit[2])) |> round()
  
  
  counts_group1 <- draw_uniform_pop(template_group1_mu, template_sigma, template_Q, libsizes_1,
                                    template_tetra_corr, template_rel_corr)
  
  counts_group2 <- draw_uniform_pop(template_group2_mu, template_sigma, template_Q, libsizes_2,
                                    template_tetra_corr, template_rel_corr)
  
  count_mat <- matrix(0, nrow=num_sample, ncol=num_taxa)
  count_mat[metadata_df$main_vec == 0, ] <- counts_group1
  count_mat[metadata_df$main_vec == 1, ] <- counts_group2
  
  return(count_mat)
  
}


