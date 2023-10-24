rm(list=ls())
suppressMessages(library(ALDEx2))

simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
aldex_folder <- file.path(simulation_folder, 'Aldex2')

source(file.path(aldex_folder, "aldex2_utils.R"))


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

prevalences <- rowMeans(simulated_data$count_mat > 0)
count_mat_subset <- simulated_data$count_mat[prevalences > 0.1, ]

ptm <- proc.time()
aldex_DA <- aldex(reads=count_mat_subset,
                                   conditions=simulated_data$sample_metadata$Group,
                                   test="t")
duration <- proc.time() - ptm


if(DAprop == 0){
  performance_aldex_DA <- evaluation(taxa_info_DA, aldex_DA, nullcase=T)
} else{
  performance_aldex_DA <- evaluation(taxa_info_DA, aldex_DA, nullcase=F)
}

performance_aldex <- evaluation(taxa_info_DA, aldex_DA)

output_aldex <- list(performance_aldex=performance_aldex,
                     duration=duration)

output_filename <- sprintf("Aldex_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
saveRDS(output_aldex, file=file.path(aldex_folder, mode, sprintf("seed_%d", sim_seed), output_filename))

