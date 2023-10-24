# library(ggplot2)
remove(list=ls())

folder <- '/home/wangmk/MDAWG/POLDA_example/simulation'
source(file.path(folder, "simulate_DirMultinom/utils.R"))

# estimate the dirichlet distribution parameter and the sequencing depth parameters
if(!file.exists(file.path(folder, "data/AGP_dirgamma.RData"))){
  AGP_template <- readRDS(file.path(folder, 'data/AGP_template.rds'))
  AGP_count_mat <- otu_table(AGP_template)@.Data
  AGP_parameters <- estimate_param(AGP_count_mat, seed=5)
  AGP_dirgamma <- AGP_parameters$dirmult_composition
  # the estimated gamma generates too much overdispersion, so I multiply all the parameters by the same constant
  scaled_AGP_dirgamma <- AGP_dirgamma * 5
  mean_depth <- AGP_parameters$Seqdepth_mean
  theta_depth <- AGP_parameters$Seqdepth_theta
  save(scaled_AGP_dirgamma, mean_depth, theta_depth, 
       file=file.path(folder, "data/AGP_dirgamma.RData"))
} else{
  load(file.path(folder, "data/AGP_dirgamma.RData"))
}


ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# null sample
# AGP_null <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
#                             seqdepth_mean = mean_depth,
#                             seqdepth_theta = theta_depth,
#                             nSample=100,
#                             DA_proportion=0, seed=ID, zinf=0.6)
# output_filename <- sprintf("AGP_simulated_%d.rds", ID)
# saveRDS(AGP_null, file=file.path(folder, 'data', 'null', output_filename))


example_data <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma/5,
              seqdepth_mean = mean_depth,
              seqdepth_theta = theta_depth,
              nSample = 100,
              grp_ratio=1,
              DA_avg_logfold = log(3),
              DA_sd_logfold=0,
              DA_proportion=0, seed=4, zinf=0,
              DA_direction="enrich",
              DA_mode="rare")



# simulations with DAs
sample_sizes <- c(50, 100)
DA_proportions <- c(5, 10, 20)
DA_directions <- c("enrich", "deplete", "mix")
DA_modes <- c("abundant", "rare")

library(tidyr)
allsettings <- crossing(sample_sizes, DA_proportions, DA_directions, DA_modes)

for (j in 1:nrow(allsettings)){
  
  samplesize <- allsettings$sample_sizes[j]
  DAprop <- allsettings$DA_proportions[j]
  DAdir <- allsettings$DA_directions[j]
  DAmode <- allsettings$DA_modes[j]
  
  AGP_simulated_data <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                      seqdepth_mean = mean_depth,
                                      seqdepth_theta = theta_depth,
                                      nSample = samplesize,
                                      grp_ratio=1,
                                      DA_avg_logfold = log(5),
                                      DA_sd_logfold=0,
                                      DA_proportion=DAprop/100, seed=ID, zinf=0.6,
                                      DA_direction=DAdir,
                                      DA_mode=DAmode)
  
  output_filename <- sprintf("AGP_simulated_%d.rds", ID)
  output_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
  saveRDS(AGP_simulated_data, file=file.path(folder, "data", output_folder, output_filename))

}


