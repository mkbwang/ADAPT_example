suppressMessages(library(pROC))

# evaluate performance
evaluation_rionorm2 <- function(taxa_truth, rionorm2_result, nullcase=FALSE){
  result <- NULL
  adj_pval <- rionorm2_result$padj
  raw_pval <- rionorm2_result$pvalue
  if (nullcase){
    FPR <- mean(raw_pval < 0.05, na.rm=T)
    result <- list(FPR = FPR, FWE=any(raw_pval < 0.05/nrow(rionorm2_result)))
  } else{
    
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    rionorm2_DA_taxa <- rionorm2_result$id[adj_pval < 0.05]
    DA_truth_binary <- taxa_truth[rionorm2_result$id, 'log_main_eff'] == 0
    roc_obj <- roc(DA_truth_binary, raw_pval)
    check_DA <- rionorm2_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DA)
    Power <- sum(check_DA) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}

