rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(dacomp))

simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
dacomp_folder <- file.path(simulation_folder, 'DACOMP')


source(file.path(dacomp_folder, "DACOMP_utils.R"))

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

prevalences <- rowMeans(simulated_data$count_mat > 0)
filter <- prevalences > 0.1
simulated_data$count_mat <- simulated_data$count_mat[filter, ]

taxa_DAinfo <- simulated_data$taxa_info

ptm <- proc.time()
dacomp_result <- pipeline(simulated_data)
duration <- proc.time() - ptm

performance_dacomp <- evaluation(taxa_DAinfo, dacomp_result, simulated_data)

output_dacomp <- list(performance_dacomp=performance_dacomp,
               duration=duration)


output_filename <- sprintf("dacomp_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
saveRDS(output_dacomp, file=file.path(dacomp_folder, mode, sprintf("seed_%d", sim_seed), output_filename))


