rm(list=ls())
library(PTDA)

simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
ptda_folder <- file.path(simulation_folder, 'PTDA')

source(file.path(ptda_folder, 'ptda_utils.R'))

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



input_filename <- sprintf("sim_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)

simulated_data <- readRDS(file.path(data_folder, mode, sprintf("seed_%d", sim_seed), input_filename))
taxa_info_DA <- simulated_data$taxa_info

prevalences <- rowMeans(simulated_data$count_mat > 0)

ptm <- proc.time()
ptda_DA_noboot <- ptda(otu_table=simulated_data$count_mat,
                  metadata=simulated_data$sample_metadata,
                  prevalence_cutoff=0.1,
                  covar="Group", genes_are_rows = TRUE, boot=F, alpha=0.05)
noboot_time <- proc.time() - ptm


ptm <- proc.time()
ptda_DA_boot <- ptda(otu_table=simulated_data$count_mat,
                       metadata=simulated_data$sample_metadata,
                       prevalence_cutoff=0.1,
                       covar="Group", genes_are_rows = TRUE, boot=T, boot_replicate = 1000, n_boot_gene = 100,
                     alpha=0.05)
boot_time <- proc.time() - ptm


true_DA_taxa <- rownames(taxa_info_DA)[taxa_info_DA$Logfold != 0]


if(DAprop == 0){
  performance_ptda_DA_noboot <- evaluation(taxa_info_DA, ptda_DA_noboot, nullcase=T)
} else{
  performance_ptda_DA_noboot <- evaluation(taxa_info_DA, ptda_DA_noboot, nullcase=F)
}


if(DAprop == 0){
  performance_ptda_DA_boot <- evaluation(taxa_info_DA, ptda_DA_boot, nullcase=T)
} else{
  performance_ptda_DA_boot <- evaluation(taxa_info_DA, ptda_DA_boot, nullcase=F)
}



output_noboot <- list(performance=performance_ptda_DA_noboot,
               ptda_result=ptda_DA_noboot,
               duration = noboot_time)

output_boot <- list(performance=performance_ptda_DA_boot,
                      ptda_result=ptda_DA_boot,
                    duration=boot_time)




output_noboot_filename <- sprintf("ptda_noboot_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
saveRDS(output_noboot, file=file.path(ptda_folder, mode, sprintf("seed_%d", sim_seed), output_noboot_filename))


output_boot_filename <- sprintf("ptda_boot_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
saveRDS(output_boot, file=file.path(ptda_folder, mode, sprintf("seed_%d", sim_seed), output_boot_filename))



