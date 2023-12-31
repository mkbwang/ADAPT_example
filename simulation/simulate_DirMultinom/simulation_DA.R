remove(list=ls())

library(dirmult)
library(optparse)
option_list = list(
  make_option(c("-n", "--nsample"), type="integer", default=100, 
              help="Number of Samples [default=100]"),
  make_option(c("-p", "--daprop"), type="double", default=0.1, 
              help="DA taxa proportion [default=0.1]"),
  make_option(c("-f", "--foldchange"), type="double", default=3, 
              help="DA taxa proportion [default=3]"),
  make_option(c("-d", "--depth"), type="integer", default=1e5,
              help="average sequencing depths [default=100000]"),
  make_option(c("-s", "--seed"), type="integer", default=1, 
              help="seed [default=1]"),
  make_option(c("-m", "--mode"), type="character", default="mix", 
              help="mode [default=mix]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

seqdepth <-opt$depth
logfold <- log(opt$foldchange)
nindv <- opt$nsample
sim_seed <- opt$seed
mode <- opt$mode
stopifnot(mode %in% c("mix", "enrich", "deplete"))
DAprop <- opt$daprop
# 
# seqdepth <-1e5
# logfold <- log(3)
# nindv <- 100
# sim_seed <- 2
# mode <- "mix"
# stopifnot(mode %in% c("mix", "enrich", "deplete"))
# DAprop <- 0.1


set.seed(sim_seed)

folder <- "/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation/simulate_DirMultinom"
source(file.path(folder, "utils.R"))

# if (!file.exists(file.path(folder, 'template_mat.rds'))){
#   adult_data <-  readRDS("/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data/global_gut/global_gut_adult.rds")
#   
#   otu_mat <- adult_data$otu_table
#   metadata <- adult_data$metadata
#   ## filter USA samples
#   US_filter <- metadata$country == "USA"
#   US_adult_otu_mat <- otu_mat[, US_filter]
#   US_metadata <- metadata[US_filter, ]
#   
#   ## retain taxa whose prevalence range between 20% to 80%
#   prevalences <- rowMeans(US_adult_otu_mat != 0)
#   taxa_filter <- prevalences > 0.2 & prevalences < 0.8
#   subset_otu_mat <- US_adult_otu_mat[taxa_filter, ]
#   
#   ## randomly select 1000 taxa
#   selected_taxa <- sample(sum(taxa_filter), size=num_taxa)
#   selected_otu_mat <- subset_otu_mat[selected_taxa, ]
#   prevalences_subset <- rowMeans(selected_otu_mat != 0)
#   saveRDS(selected_otu_mat, file=file.path(folder, 'template_mat.rds'))
# } else{
#   selected_otu_mat <- readRDS(file.path(folder, 'template_mat.rds'))
# }
# 
# 
# if(!file.exists(file.path(folder, 'template_dirparam.rds'))){
#   dirichlet_param <- dirmult(t(selected_otu_mat))
#   gamma_param <- dirichlet_param$gamma
#   saveRDS(gamma_param, file.path(folder, "template_dirparam.rds"))
# } else{
#   gamma_param <- readRDS(file.path(folder, "template_dirparam.rds"))
# }
# 
# 
# if (!file.exists(file.path(folder, 'AGP_template_dirparam.rds'))){
#   count_mat <- readRDS("/scratch/ligen_root/ligen0/wangmk/POLDA_example/real_data/AGP/AGP_template.rds")
#   dirmult_fit <- dirmult(data=t(count_mat))
#   gamma_param <- dirmult_fit$gamma
#   saveRDS(gamma_param, file.path(output_folder, "simulation/simulate_DirMultinom/AGP_template_dirparam.rds"))
# } else{
#   gamma_param <- readRDS(file.path(folder, "AGP_template_dirparam.rds"))
# }

gamma_param <- readRDS(file.path(folder, "global_gut_template_dirparam.rds"))

# simdata <- SimCount2(dirmult_composition = gamma_param, seqdepth_mean=seqdepth,
#                          seqdepth_theta=10, nSample=nindv, DA_proportion = DAprop, mean_logfold = logfold,
#                          sd_logfold=0,
#                          DA_direction = mode,
#                          seed=sim_seed)




simdata <- SimulateCount(dirmult_composition = gamma_param, seqdepth_mean=seqdepth,
                         seqdepth_theta=10, nSample=nindv, DA_proportion = DAprop, mean_logfold = logfold,
                         sd_logfold=0,
                         DA_direction = mode,
                         seed=sim_seed)



output_file <- sprintf("sim_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
saveRDS(simdata, file.path(folder, mode, sprintf("seed_%d", sim_seed), output_file))

