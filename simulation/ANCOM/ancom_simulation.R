rm(list=ls())
suppressMessages(library(ANCOMBC))
suppressMessages(library(phyloseq))


simulation_folder <- '/scratch/ligen_root/ligen0/wangmk/PTDA_example/simulation'
data_folder <- file.path(simulation_folder, 'simulate_DirMultinom')
ancom_folder <- file.path(simulation_folder, 'ANCOM')


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
out <- ancom(data = phyobj,
             p_adj_method = "BH", prv_cut = 0.1, lib_cut = 0,
             main_var = "Group", struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05, n_cl = 4)
duration <- proc.time() - ptm

result <- out$res
true_DA_taxa <- rownames(taxa_info_DA)[taxa_info_DA$Logfold != 0]
ancom_DA_taxa <- result$taxon[result$detected_0.7]
check_DAtaxa <- ancom_DA_taxa %in% true_DA_taxa
FDR <- 1 - mean(check_DAtaxa)
Power <- sum(check_DAtaxa) / length(true_DA_taxa)
performance_ancom <- list(FDR = FDR,
               Power=Power)

output_ancom <- list(performance_ancom=performance_ancom,
                     duration=duration)

output_filename <- sprintf("ANCOM_n%d_p%d_d%d_f%.1f_s%d.rds", nindv, DAprop*100, seqdepth, exp(logfold), sim_seed)
saveRDS(output_ancom, file=file.path(ancom_folder, mode, sprintf("seed_%d", sim_seed), output_filename))


