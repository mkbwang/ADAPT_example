suppressMessages(library(pROC))

# evaluate performance
evaluation_raida <- function(taxa_truth, raida_result, nullcase=FALSE){
  result <- NULL
  adj_pval <- raida_result$p.adj
  raw_pval <- raida_result$p
  if (nullcase){
    FPR <- mean(raw_pval < 0.05, na.rm=T)
    result <- list(FPR = FPR, FWE=any(raw_pval < 0.05/nrow(raida_result)))
  } else{
    
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    raida_DA_taxa <- rownames(raida_result)[adj_pval < 0.05]
    DA_truth_binary <- taxa_truth[rownames(raida_result), 'log_main_eff'] == 0
    roc_obj <- suppressMessages(roc(DA_truth_binary, raw_pval))
    check_DA <- raida_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DA)
    Power <- sum(check_DA) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}

