remove(list=ls())

folder <- '/home/wangmk/UM/Research/MDAWG/POLDA_example/simulation'
source(file.path(folder, "simulate_DirMultinom/utils.R"))


rm(list=ls())
adult_data <-  readRDS("real_data/Yatsunenko_data/global_gut_adult.rds")

otu_mat <- adult_data$otu_table
metadata <- adult_data$metadata
## filter USA samples
US_filter <- metadata$country == "USA"
US_adult_otu_mat <- otu_mat[, US_filter]
US_metadata <- metadata[US_filter, ]

## retain taxa whose prevalence range between 20% to 80%
prevalences <- rowMeans(US_adult_otu_mat != 0)
taxa_filter <- prevalences > 0.2 & prevalences < 0.8
subset_otu_mat <- US_adult_otu_mat[taxa_filter, ]



set.seed(1)
num_taxa <- 1000
## randomly select 2000 taxa
selected_taxa <- sample(sum(taxa_filter), size=num_taxa)
selected_otu_mat <- subset_otu_mat[selected_taxa, ]
## estimate Dirichlet distribution parameters
dirichlet_param <- dirmult(t(selected_otu_mat))
gamma_param <- sort(dirichlet_param$gamma)

## balanced situation (swap gamma parameters between group1 and group2)

gen_dirparam <- function(baseline, propDA=0.1, balanced=T, enrich=F, seed=1){
  stopifnot(propDA < 0.2)
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
    } else{
      g2_gamma[DA_prevalenttaxa_ID] <- sample(raretaxa_gamma, size=num_DAtaxa)
    }
  }
  result <- list(g1_gamma=g1_gamma, g2_gamma=g2_gamma, log_fold_change = log(g2_gamma/g1_gamma))
  return(result)
}


trial_balanced <- gen_dirparam(baseline=gamma_param, propDA=0.1, balanced=T, seed=5)
trial_enriched <- gen_dirparam(baseline=gamma_param, propDA=0.1, balanced=F, enrich = T, seed=8)
trial_depleted <- gen_dirparam(baseline=gamma_param, propDA=0.1, balanced=F, enrich = F, seed=6)

## given number of samples and average sequencing depths, generate count matrices

