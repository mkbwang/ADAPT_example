
rm(list=ls())
LOCOM_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/LOCOM"

sample_sizes <- c(50, 100)
DA_proportions <- c(5, 10, 20)
DA_directions <- c("enrich", "deplete", "mix")
DA_modes <- c("rare", "abundant")

library(tidyr)
allsettings <- crossing(sample_sizes, DA_proportions, DA_directions, DA_modes)



locom_performances_allsettings <- list()

for (j in 1:nrow(allsettings)){
  samplesize <- allsettings$sample_sizes[j]
  DAprop <- allsettings$DA_proportions[j]
  DAdir <- allsettings$DA_directions[j]
  DAmode <- allsettings$DA_modes[j]
  
  input_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
  allfiles <- list.files(file.path(LOCOM_folder, input_folder)) |> 
    tools::file_path_sans_ext()
  available_IDs <- sapply(allfiles, function(name) strsplit(name, split="_")[[1]][3]) |>
    as.integer()
  
  locom_performance <- data.frame(ID=available_IDs,
                                  reftaxa_error=rep(0, length(available_IDs)),
                                  AUROC=rep(0, length(available_IDs)),
                                  FDR=rep(0, length(available_IDs)),
                                  Power=rep(0, length(available_IDs)))
  locom_performance$sample_size <- samplesize
  locom_performance$Proportion <- DAprop
  locom_performance$Direction <- DAdir
  locom_performance$Type <- DAmode
  
  
  for (k in 1:length(available_IDs)){
    ID <- available_IDs[k]
    input_filename <- sprintf("locom_summary_%d.rds", ID)
    model_summary <- readRDS(file.path(LOCOM_folder, input_folder, input_filename))
    model_performance <- model_summary$performance
    locom_performance$Power[k] <- model_performance$Power
    # locom_performance$reftaxa_error[ID] <- model_performance$reftaxa_error
    locom_performance$AUROC[k] <- model_performance$AUROC
    locom_performance$FDR[k] <- model_performance$FDR
  }
  locom_performances_allsettings[[j]] <- locom_performance
}

locom_new_results <- do.call(rbind, locom_performances_allsettings)
locom_new_results$Method <- "LOCOM"
saveRDS(locom_new_results, file=file.path(LOCOM_folder, "locom_performance.rds"))
