suppressMessages(library(pROC))

# evaluate performance
evaluation_rdb <- function(taxa_truth, rdb_result, nullcase=FALSE){
  result <- NULL
  if (nullcase){
    result <- list(FWE = any(rdb_result, na.rm=T))
  } else{
    
    DA_truth_binary <- taxa_truth$log_main_eff != 0
    check_DAtaxa <- rdb_result & DA_truth_binary
    FDR <- 1 - sum(check_DAtaxa)/sum(rdb_result)
    Power <- sum(check_DAtaxa) / sum(DA_truth_binary)
    result <- list(FDR = FDR,
                   FPR=NA,
                   Power=Power)
  }
  
  return(result)
}

