suppressMessages(library(pROC))



# evaluate performance
evaluation_ancom <- function(taxa_truth, ancom_result, nullcase=FALSE){
  result <- NULL
  if (nullcase){
    result <- list(FPR=NA,
                   FWE=NA)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    ancom_W <- ancom_result$W
    roc_obj <-  suppressMessages(roc(!DA_truth_binary, ancom_W))
    ancom_DA_taxa <- ancom_result$taxon[ancom_result$detected_0.7]
    check_DAtaxa <- ancom_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   FPR=NA,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


