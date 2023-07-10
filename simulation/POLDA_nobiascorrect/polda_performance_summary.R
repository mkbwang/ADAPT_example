

POLDA_folder <- "/home/wangmk/MDAWG/POLDA_example/simulation/POLDA_nobiascorrect"

sample_sizes <- c(50, 100)
DA_proportions <- c(5, 10, 20)
DA_directions <- c("enrich", "deplete", "mix")
DA_modes <- c("rare", "abundant")

library(tidyr)
allsettings <- crossing(sample_sizes, DA_proportions, DA_directions, DA_modes)

IDs <- setdiff(seq(1, 100),
               c(3,14,18,21,24,36,43,
                 45,52,56,57,61,66,71,73,
                 83,87,89,91,94,95,98)) # some files are old

polda_performances_allsettings <- list()


for (j in 1:nrow(allsettings)){
  samplesize <- allsettings$sample_sizes[j]
  DAprop <- allsettings$DA_proportions[j]
  DAdir <- allsettings$DA_directions[j]
  DAmode <- allsettings$DA_modes[j]
  
  input_folder <- sprintf("N%d_P%d/%s_%s", samplesize, DAprop, DAdir, DAmode)
  allfiles <- list.files(file.path(POLDA_folder, input_folder)) |> 
    tools::file_path_sans_ext()
  available_IDs <- sapply(allfiles, function(name) strsplit(name, split="_")[[1]][3]) |>
    as.integer()
  usable_IDs <- intersect(IDs, available_IDs)
  
  polda_performance <- data.frame(ID=usable_IDs,
                                  reftaxa_error=rep(0, length(usable_IDs)),
                                  AUROC=rep(0, length(usable_IDs)),
                                  FDR=rep(0, length(usable_IDs)),
                                  Power=rep(0, length(usable_IDs)))
  
  polda_performance$sample_size <- samplesize
  polda_performance$Proportion <- DAprop
  polda_performance$Direction <- DAdir
  polda_performance$Type <- DAmode
  
  for (k in 1:length(usable_IDs)){
    ID <- usable_IDs[k]
    input_filename <- sprintf("polda_summary_%d.rds", ID)
    model_summary <- readRDS(file.path(POLDA_folder, input_folder, input_filename))
    model_performance <- model_summary$performance
    polda_performance$Power[k] <- model_performance$Power
    polda_performance$reftaxa_error[k] <- model_performance$reftaxa_error
    polda_performance$AUROC[k] <- model_performance$AUROC
    polda_performance$FDR[k] <- model_performance$FDR
  }
  polda_performances_allsettings[[j]] <- polda_performance
}

polda_old_results <- do.call(rbind, polda_performances_allsettings)
polda_old_results$Method <- "POLDA_Old"
saveRDS(polda_old_results, file=file.path(POLDA_folder, "polda_old_performance.rds"))
