rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(LOCOM))

simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
locom_folder <- file.path(simulation_folder, 'LOCOM')

source(file.path(locom_folder, "locom_utils.R"))


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


# seqdepth <-1e5
# nindv <- 50
# sim_seed <- 30
# mode <- "enrich"
# DAprop <- 0.1


input_filename <- sprintf("sim_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
simulated_data <- readRDS(file.path(data_folder, mode, sprintf("seed_%d", sim_seed), input_filename))

taxa_info_DA <- simulated_data$taxa_info

ptm <- proc.time()
locom_result <- locom(otu.table=t(simulated_data$count_mat),
                  Y=simulated_data$sample_metadata$Group, fdr.nominal=0.05, prev.cut=0.1,
                  n.rej.stop=100, seed=1, n.perm.max=2e4, Firth.thresh = 1,
                  n.cores=4)
duration <- proc.time() - ptm

performance_locom <- evaluation(taxa_info_DA, locom_result, nullcase=F)

output_locom <- list(performance_locom=performance_locom,
                     duration=duration)

output_filename <- sprintf("LOCOM_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), 
                           sim_seed)
saveRDS(output_locom, file=file.path(locom_folder, mode, sprintf("seed_%d", sim_seed), output_filename))


