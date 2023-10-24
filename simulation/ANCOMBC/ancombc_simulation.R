rm(list=ls())
suppressMessages(library(ANCOMBC))
suppressMessages(library(phyloseq))


simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
ancombc_folder <- file.path(simulation_folder, 'ANCOMBC')

source(file.path(ancombc_folder, "ancombc_utils.R"))


library(optparse)
option_list = list(
  make_option(c("-n", "--nsample"), type="integer", default=100, 
              help="Number of Samples [default=100]"),
  make_option(c("-p", "--daprop"), type="double", default=0.1, 
              help="DA taxa proportion [default=0.1]"),
  make_option(c("-d", "--depth"), type="integer", default=1e5,
              help="average sequencing depths [default=100,000]"),
  make_option(c("-s", "--seed"), type="integer", default=1, 
              help="seed [default=1]"),
  make_option(c("-f", "--foldchange"), type="double", default=3, 
              help="DA taxa proportion [default=3]"),
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


# seqdepth <-2e4
# nindv <- 100
# sim_seed <- 3
# mode <- "mix"
# DAprop <- 0.1


input_filename <- sprintf("sim_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
simulated_data <- readRDS(file.path(data_folder, mode, sprintf("seed_%d", sim_seed), input_filename))

taxa_info_DA <- simulated_data$taxa_info

phyobj <- phyloseq(otu_table(simulated_data$count_mat, taxa_are_rows = TRUE),
                   sample_data(simulated_data$sample_metadata))

ptm <- proc.time()
ancombc_result <- ancombc(phyobj, formula = 'Group', p_adj_method='BH',
                                        group='Group', struc_zero=TRUE, prv_cut=0.1, n_cl=3)
duration <- proc.time() - ptm

performance_ancombc <- evaluation(taxa_info_DA,
                                               ancombc_result)

output_ancombc <- list(performance_ancombc=performance_ancombc,
                       duration=duration)

output_filename <- sprintf("ANCOMBC_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
saveRDS(output_ancombc, file=file.path(ancombc_folder, mode, sprintf("seed_%d", sim_seed), output_filename))


