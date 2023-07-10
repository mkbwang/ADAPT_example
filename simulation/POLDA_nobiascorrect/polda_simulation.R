rm(list=ls())
library(POLDA)

args <- commandArgs(trailingOnly = TRUE)
bias_correction <- as.logical(args[1]) # whether to use bias correction

data_folder <- '/home/wangmk/MDAWG/POLDA_example/simulation/data'
if (bias_correction){
  polda_folder <- '/home/wangmk/MDAWG/POLDA_example/simulation/POLDA'
} else{
  polda_folder <- '/home/wangmk/MDAWG/POLDA_example/simulation/POLDA_nobiascorrect'
}

source(file.path(polda_folder, 'polda_utils.R'))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
polda_null <- polda(otu_table=AGP_null$count_mat,
                    metadata=AGP_null$sample_metadata,
                    prevalence_cutoff=0.05,
                    covar="X",
                    bias_correct = bias_correction)
performance_polda_null <- evaluation(taxa_truth = taxa_info_null,
                                     polda_result = polda_null,
                                     nullcase=TRUE)
output_null <- list(performance=performance_polda_null,
                                   polda_result=polda_null)
saveRDS(output_null,
        file.path(polda_folder, "null",
                  sprintf("summary_null_%d.rds", ID)))



# simulations with DAs
sample_sizes <- c(50, 100)
DA_proportions <- c(5, 10, 20)
DA_directions <- c("enrich", "deplete", "mix")
DA_modes <- c("rare", "abundant")

library(tidyr)
allsettings <- crossing(sample_sizes, DA_proportions, DA_directions, DA_modes)

for (j in 1:nrow(allsettings)){
  
  samplesize <- allsettings$sample_sizes[j]
  DAprop <- allsettings$DA_proportions[j]
  DAdir <- allsettings$DA_directions[j]
  DAmode <- allsettings$DA_modes[j]
  
  input_filename <- sprintf("AGP_simulated_%d.rds", ID)
  input_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
  AGP_simulated_data <- readRDS(file.path(data_folder, input_folder, input_filename))
  
  taxa_info_DA <- AGP_simulated_data$taxa_info
  polda_DA <- polda(otu_table=AGP_simulated_data$count_mat,
                    metadata=AGP_simulated_data$sample_metadata,
                    prevalence_cutoff=0.05,
                    covar="X",
                    bias_correct = bias_correction)
  performance_polda_DA <- evaluation(taxa_truth = taxa_info_DA,
                                     polda_result = polda_DA,
                                     nullcase=FALSE)
  
  output <- list(performance=performance_polda_DA,
                 polda_result=polda_DA)
  
  output_filename <- sprintf("polda_summary_%d.rds", ID)
  output_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
  saveRDS(output, file=file.path(polda_folder, output_folder, output_filename))
  
}
  



# 
# 
# # rare DA taxa, low proportion(5%)
# AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
#                                              sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
# taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
# polda_unbalanced_rare_low <- polda(otu_table = AGP_unbalanced_rare_low$count_mat,
#                                    metadata = AGP_unbalanced_rare_low$sample_metadata,
#                                    covar="X")
# performance_unbalanced_rare_low <- evaluation(taxa_truth = taxa_info_unbalanced_rare_low,
#                                               polda_result = polda_unbalanced_rare_low)
# output_unbalanced_rare_low <- list(performance=performance_unbalanced_rare_low,
#                                    polda_result=polda_unbalanced_rare_low)
# saveRDS(output_unbalanced_rare_low, 
#         file.path(polda_folder, "DA_rare", "low", 
#                   sprintf("summary_unbalanced_rare_low_%d.rds", ID)))
# 
# 
# # rare DA taxa, medium proportion(10%)
# AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
#                                              sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
# taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
# polda_unbalanced_rare_medium <- polda(otu_table = AGP_unbalanced_rare_medium$count_mat,
#                                    metadata = AGP_unbalanced_rare_medium$sample_metadata,
#                                    covar="X")
# performance_unbalanced_rare_medium <- evaluation(taxa_truth = taxa_info_unbalanced_rare_medium,
#                                               polda_result = polda_unbalanced_rare_medium)
# output_unbalanced_rare_medium <- list(performance=performance_unbalanced_rare_medium,
#                                    polda_result=polda_unbalanced_rare_medium)
# saveRDS(output_unbalanced_rare_medium, 
#         file.path(polda_folder, "DA_rare", "medium", 
#                   sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))
# 
# 
# # rare DA taxa, high proportion(20%)
# AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
#                                                 sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
# taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
# polda_unbalanced_rare_high <- polda(otu_table = AGP_unbalanced_rare_high$count_mat,
#                                       metadata = AGP_unbalanced_rare_high$sample_metadata,
#                                       covar="X")
# performance_unbalanced_rare_high <- evaluation(taxa_truth = taxa_info_unbalanced_rare_high,
#                                                  polda_result = polda_unbalanced_rare_high)
# output_unbalanced_rare_high <- list(performance=performance_unbalanced_rare_high,
#                                       polda_result=polda_unbalanced_rare_high)
# saveRDS(output_unbalanced_rare_high, 
#         file.path(polda_folder, "DA_rare", "high", 
#                   sprintf("summary_unbalanced_rare_high_%d.rds", ID)))
# 
# 
# 
# # abundant DA taxa, low proportion(5%)
# AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
#                                               sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
# taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
# polda_unbalanced_abundant_low <- polda(otu_table = AGP_unbalanced_abundant_low$count_mat,
#                                     metadata = AGP_unbalanced_abundant_low$sample_metadata,
#                                     covar="X")
# performance_unbalanced_abundant_low <- evaluation(taxa_truth = taxa_info_unbalanced_abundant_low,
#                                                polda_result = polda_unbalanced_abundant_low)
# output_unbalanced_abundant_low <- list(performance=performance_unbalanced_abundant_low,
#                                     polda_result=polda_unbalanced_abundant_low)
# saveRDS(output_unbalanced_abundant_low, 
#         file.path(polda_folder, "DA_abundant", "low", 
#                   sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))
# 
# 
# 
# # abundant DA taxa, medium proportion(10%)
# AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
#                                                  sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
# taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
# polda_unbalanced_abundant_medium <- polda(otu_table = AGP_unbalanced_abundant_medium$count_mat,
#                                        metadata = AGP_unbalanced_abundant_medium$sample_metadata,
#                                        covar="X")
# performance_unbalanced_abundant_medium <- evaluation(taxa_truth = taxa_info_unbalanced_abundant_medium,
#                                                   polda_result = polda_unbalanced_abundant_medium)
# output_unbalanced_abundant_medium <- list(performance=performance_unbalanced_abundant_medium,
#                                        polda_result=polda_unbalanced_abundant_medium)
# saveRDS(output_unbalanced_abundant_medium, 
#         file.path(polda_folder, "DA_abundant", "medium", 
#                   sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))
# 
# 
# # abundant DA taxa, high proportion(20%)
# AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
#                                                     sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
# taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
# polda_unbalanced_abundant_high <- polda(otu_table = AGP_unbalanced_abundant_high$count_mat,
#                                           metadata = AGP_unbalanced_abundant_high$sample_metadata,
#                                           covar="X")
# performance_unbalanced_abundant_high <- evaluation(taxa_truth = taxa_info_unbalanced_abundant_high,
#                                                      polda_result = polda_unbalanced_abundant_high)
# output_unbalanced_abundant_high <- list(performance=performance_unbalanced_abundant_high,
#                                           polda_result=polda_unbalanced_abundant_high)
# saveRDS(output_unbalanced_abundant_high, 
#         file.path(polda_folder, "DA_abundant", "high", 
#                   sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))
# 


# mixed DA taxa, low proportion(5%)
# AGP_unbalanced_mixed_low <- readRDS(file.path(data_folder, "DA_mixed", "low",
#                                               sprintf("AGP_simulation_unbalanced_mixed_low_%d.rds", ID)))
# taxa_info_unbalanced_mixed_low <- AGP_unbalanced_mixed_low$taxa_info
# polda_unbalanced_mixed_low <- polda(otu_table = AGP_unbalanced_mixed_low$count_mat,
#                                     metadata = AGP_unbalanced_mixed_low$sample_metadata,
#                                     covar="X")
# performance_unbalanced_mixed_low <- evaluation(taxa_truth = taxa_info_unbalanced_mixed_low,
#                                                polda_result = polda_unbalanced_mixed_low)
# output_unbalanced_mixed_low <- list(performance=performance_unbalanced_mixed_low,
#                                     polda_result=polda_unbalanced_mixed_low)
# saveRDS(output_unbalanced_mixed_low, 
#         file.path(polda_folder, "DA_mixed", "low", 
#                   sprintf("summary_unbalanced_mixed_low_%d.rds", ID)))
# 
# 
# 
# # mixed DA taxa, medium proportion(10%)
# AGP_unbalanced_mixed_medium <- readRDS(file.path(data_folder, "DA_mixed", "medium",
#                                                  sprintf("AGP_simulation_unbalanced_mixed_medium_%d.rds", ID)))
# taxa_info_unbalanced_mixed_medium <- AGP_unbalanced_mixed_medium$taxa_info
# polda_unbalanced_mixed_medium <- polda(otu_table = AGP_unbalanced_mixed_medium$count_mat,
#                                        metadata = AGP_unbalanced_mixed_medium$sample_metadata,
#                                        covar="X")
# performance_unbalanced_mixed_medium <- evaluation(taxa_truth = taxa_info_unbalanced_mixed_medium,
#                                                   polda_result = polda_unbalanced_mixed_medium)
# output_unbalanced_mixed_medium <- list(performance=performance_unbalanced_mixed_medium,
#                                        polda_result=polda_unbalanced_mixed_medium)
# saveRDS(output_unbalanced_mixed_medium, 
#         file.path(polda_folder, "DA_mixed", "medium", 
#                   sprintf("summary_unbalanced_mixed_medium_%d.rds", ID)))
# 
# 
# # mixed DA taxa, high proportion(20%)
# AGP_unbalanced_mixed_high <- readRDS(file.path(data_folder, "DA_mixed", "high",
#                                                sprintf("AGP_simulation_unbalanced_mixed_high_%d.rds", ID)))
# taxa_info_unbalanced_mixed_high <- AGP_unbalanced_mixed_high$taxa_info
# polda_unbalanced_mixed_high <- polda(otu_table = AGP_unbalanced_mixed_high$count_mat,
#                                      metadata = AGP_unbalanced_mixed_high$sample_metadata,
#                                      covar="X")
# performance_unbalanced_mixed_high <- evaluation(taxa_truth = taxa_info_unbalanced_mixed_high,
#                                                 polda_result = polda_unbalanced_mixed_high)
# output_unbalanced_mixed_high <- list(performance=performance_unbalanced_mixed_high,
#                                      polda_result=polda_unbalanced_mixed_high)
# saveRDS(output_unbalanced_mixed_high, 
#         file.path(polda_folder, "DA_mixed", "high", 
#                   sprintf("summary_unbalanced_mixed_high_%d.rds", ID)))




