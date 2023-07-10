

POLDA_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/POLDA"

sample_sizes <- c(50, 100)
DA_proportions <- c(5, 10, 20)
DA_directions <- c("enrich", "deplete", "mix")
DA_modes <- c("rare", "abundant")

library(tidyr)
allsettings <- crossing(sample_sizes, DA_proportions, DA_directions, DA_modes)



polda_performances_allsettings <- list()

for (j in 1:nrow(allsettings)){
  samplesize <- allsettings$sample_sizes[j]
  DAprop <- allsettings$DA_proportions[j]
  DAdir <- allsettings$DA_directions[j]
  DAmode <- allsettings$DA_modes[j]
  
  polda_performance <- data.frame(ID=seq(1, 100),
                                  reftaxa_error=rep(0, 100),
                                  AUROC=rep(0, 100),
                                  FDR=rep(0, 100),
                                  Power=rep(0, 100))
  polda_performance$sample_size <- samplesize
  polda_performance$Proportion <- DAprop
  polda_performance$Direction <- DAdir
  polda_performance$Type <- DAmode
  
  input_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
  for (ID in 1:100){
    input_filename <- sprintf("polda_summary_%d.rds", ID)
    model_summary <- readRDS(file.path(POLDA_folder, input_folder, input_filename))
    model_performance <- model_summary$performance
    polda_performance$Power[ID] <- model_performance$Power
    polda_performance$reftaxa_error[ID] <- model_performance$reftaxa_error
    polda_performance$AUROC[ID] <- model_performance$AUROC
    polda_performance$FDR[ID] <- model_performance$FDR
  }
  polda_performances_allsettings[[j]] <- polda_performance
}

polda_new_results <- do.call(rbind, polda_performances_allsettings)
saveRDS(polda_new_results, file=file.path(POLDA_folder, "polda_new_performance.rds"))
