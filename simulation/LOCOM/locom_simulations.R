rm(list=ls())
library(dplyr)
library(LOCOM)

data_folder <- '/home/wangmk/MDAWG/POLDA_example/simulation/data'
locom_folder <- '/home/wangmk/MDAWG/POLDA_example/simulation/LOCOM'


source(file.path(locom_folder, "locom_utils.R"))
ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 24


# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata <- AGP_null$sample_metadata
count_mat_null <- AGP_null$count_mat
locom_null <- locom(otu.table=t(count_mat_null),
                  Y=sample_metadata$X, fdr.nominal=0.05, prev.cut=0.05,
                  n.rej.stop=100, seed=1, n.perm.max=1e4,
                  n.cores=2)
performance_null <- evaluation(taxa_info_null, locom_null, nullcase=TRUE)
output_null <- list(performance=performance_null,
                    locom_result=locom_null)
saveRDS(output_null, 
        file.path(locom_folder, "null", 
                  sprintf("summary_null_%d.rds", ID)))


# simulations with DAs
sample_sizes <- c(50, 100)
DA_proportions <- c(5, 10, 20)
DA_directions <- c("enrich", "deplete", "mix")
DA_modes <- c("rare", "abundant") # "abundant",

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
    
    locom_DA <- locom(otu.table=t(AGP_simulated_data$count_mat),
                      Y=AGP_simulated_data$sample_metadata$X, fdr.nominal=0.05, prev.cut=0.05,
                      n.rej.stop=100, seed=1, n.perm.max=1e4,
                      n.cores=4)
    
    performance_locom_DA <- evaluation(taxa_truth = taxa_info_DA,
                                       locom_result = locom_DA,
                                       nullcase=FALSE)
    
    output_DA <- list(performance=performance_locom_DA,
                                       locom_result=locom_DA)
    
    output_filename <- sprintf("locom_summary_%d.rds", ID)
    output_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
    saveRDS(output_DA, file=file.path(locom_folder, output_folder, output_filename))
    
}









# # rare DA taxa, low proportion(5%)
# AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
#                                              sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
# taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
# count_mat_unbalanced_rare_low <- AGP_unbalanced_rare_low$count_mat
# sample_metadata <- AGP_unbalanced_rare_low$sample_metadata
# locom_unbalanced_rare_low <- locom(otu.table=t(count_mat_unbalanced_rare_low),
#                                    Y=sample_metadata$X, fdr.nominal=0.05, prev.cut=0.1,
#                                    n.rej.stop=2, seed=1,
#                                    n.cores=1)
# performance_unbalanced_rare_low <- evaluation(taxa_info_unbalanced_rare_low,
#                                               locom_unbalanced_rare_low)
# output_unbalanced_rare_low <- list(performance=performance_unbalanced_rare_low,
#                                    locom_result=locom_unbalanced_rare_low)
# saveRDS(output_unbalanced_rare_low, 
#         file.path(locom_folder, "DA_rare", "low", 
#                   sprintf("summary_unbalanced_rare_low_%d.rds", ID)))
# 
# 
# # rare DA taxa, medium proportion(10%)
# AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
#                                                 sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
# taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
# sample_metadata <- AGP_unbalanced_rare_medium$sample_metadata
# count_mat_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$count_mat
# locom_unbalanced_rare_medium <- locom(otu.table=t(count_mat_unbalanced_rare_medium),
#                                    Y=sample_metadata$X, fdr.nominal=0.05, prev.cut=0.1,
#                                    n.rej.stop=5, seed=1,
#                                    n.cores=2)
# performance_unbalanced_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium,
#                                                  locom_unbalanced_rare_medium)
# output_unbalanced_rare_medium <- list(performance=performance_unbalanced_rare_medium,
#                                    locom_result=locom_unbalanced_rare_medium)
# saveRDS(output_unbalanced_rare_medium, 
#         file.path(locom_folder, "DA_rare", "medium", 
#                   sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))
# 
# 
# # rare DA taxa, high proportion(20%)
# AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
#                                                 sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
# taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
# sample_metadata <- AGP_unbalanced_rare_high$sample_metadata
# count_mat_unbalanced_rare_high <- AGP_unbalanced_rare_high$count_mat
# locom_unbalanced_rare_high <- locom(otu.table=t(count_mat_unbalanced_rare_high),
#                                       Y=sample_metadata$X, fdr.nominal=0.05, prev.cut=0.1,
#                                       n.rej.stop=5, seed=1,
#                                       n.cores=2)
# performance_unbalanced_rare_high <- evaluation(taxa_info_unbalanced_rare_high,
#                                                  locom_unbalanced_rare_high)
# output_unbalanced_rare_high <- list(performance=performance_unbalanced_rare_high,
#                                       locom_result=locom_unbalanced_rare_high)
# saveRDS(output_unbalanced_rare_high, 
#         file.path(locom_folder, "DA_rare", "high", 
#                   sprintf("summary_unbalanced_rare_high_%d.rds", ID)))
# 
# 
# # abundant DA taxa, low proportion(5%)
# AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
#                                                  sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
# taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
# count_mat_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$count_mat
# sample_metadata <- AGP_unbalanced_abundant_low$sample_metadata
# locom_unbalanced_abundant_low <- locom(otu.table=t(count_mat_unbalanced_abundant_low),
#                                        Y=sample_metadata$X, fdr.nominal=0.05, prev.cut=0.1,
#                                        n.rej.stop=5, seed=1,
#                                        n.cores=2)
# performance_unbalanced_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low,
#                                                   locom_unbalanced_abundant_low)
# output_unbalanced_abundant_low <- list(performance=performance_unbalanced_abundant_low,
#                                        locom_result=locom_unbalanced_abundant_low)
# saveRDS(output_unbalanced_abundant_low, 
#         file.path(locom_folder, "DA_abundant", "low", 
#                   sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))
# 
# 
# # abundant DA taxa, medium proportion(10%)
# AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
#                                                     sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
# taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
# sample_metadata <- AGP_unbalanced_abundant_medium$sample_metadata
# count_mat_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$count_mat
# locom_unbalanced_abundant_medium <- locom(otu.table=t(count_mat_unbalanced_abundant_medium),
#                                           Y=sample_metadata$X, fdr.nominal=0.05, prev.cut=0.1,
#                                           n.rej.stop=5, seed=1,
#                                           n.cores=2)
# performance_unbalanced_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium,
#                                                      locom_unbalanced_abundant_medium)
# output_unbalanced_abundant_medium <- list(performance=performance_unbalanced_abundant_medium,
#                                           locom_result=locom_unbalanced_abundant_medium)
# saveRDS(output_unbalanced_abundant_medium, 
#         file.path(locom_folder, "DA_abundant", "medium", 
#                   sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))
# 
# 
# # abundant DA taxa, high proportion(20%)
# AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
#                                                   sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
# taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
# sample_metadata <- AGP_unbalanced_abundant_high$sample_metadata
# count_mat_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$count_mat
# locom_unbalanced_abundant_high <- locom(otu.table=t(count_mat_unbalanced_abundant_high),
#                                         Y=sample_metadata$X, fdr.nominal=0.05, prev.cut=0.1,
#                                         n.rej.stop=5, seed=1,
#                                         n.cores=2)
# performance_unbalanced_abundant_high <- evaluation(taxa_info_unbalanced_abundant_high,
#                                                    locom_unbalanced_abundant_high)
# output_unbalanced_abundant_high <- list(performance=performance_unbalanced_abundant_high,
#                                         locom_result=locom_unbalanced_abundant_high)
# saveRDS(output_unbalanced_abundant_high, 
#         file.path(locom_folder, "DA_abundant", "high", 
#                   sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))
