

metagenomeseq_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/MetagenomeSeq"

sample_sizes <- c(50, 100)
DA_proportions <- c(5, 10, 20)
DA_directions <- c("enrich", "deplete", "mix")
DA_modes <- c("rare", "abundant")

library(tidyr)
allsettings <- crossing(sample_sizes, DA_proportions, DA_directions, DA_modes)



metagenomeseq_performances_allsettings <- list()

for (j in 1:nrow(allsettings)){
  samplesize <- allsettings$sample_sizes[j]
  DAprop <- allsettings$DA_proportions[j]
  DAdir <- allsettings$DA_directions[j]
  DAmode <- allsettings$DA_modes[j]
  
  metagenomeseq_performance <- data.frame(ID=seq(1, 100),
                                  reftaxa_error=rep(0, 100),
                                  AUROC=rep(0, 100),
                                  FDR=rep(0, 100),
                                  Power=rep(0, 100))
  metagenomeseq_performance$sample_size <- samplesize
  metagenomeseq_performance$Proportion <- DAprop
  metagenomeseq_performance$Direction <- DAdir
  metagenomeseq_performance$Type <- DAmode
  
  input_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
  for (ID in 1:100){
    input_filename <- sprintf("metagenomeSeq_summary_%d.rds", ID)
    model_summary <- readRDS(file.path(metagenomeseq_folder, input_folder, input_filename))
    model_performance <- model_summary$performance
    metagenomeseq_performance$Power[ID] <- model_performance$Power
    # polda_performance$reftaxa_error[ID] <- model_performance$reftaxa_error
    metagenomeseq_performance$AUROC[ID] <- model_performance$AUROC
    metagenomeseq_performance$FDR[ID] <- model_performance$FDR
  }
  metagenomeseq_performances_allsettings[[j]] <- metagenomeseq_performance
}

metagenomeseq_results <- do.call(rbind, metagenomeseq_performances_allsettings)
metagenomeseq_results$Method <- "MetagenomeSeq"
saveRDS(metagenomeseq_results, file=file.path(metagenomeseq_folder, "metagenomeseq_performance.rds"))
