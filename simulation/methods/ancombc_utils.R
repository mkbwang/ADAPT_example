suppressMessages(library(pROC))



# evaluate performance
evaluation_ancombc <- function(taxa_truth, ancombc_result, nullcase=FALSE){
  result <- NULL
  output_df <- ancombc_result$res
  pvals <- output_df$p_exposure1
  qvals <- output_df$q_exposure1
  if (nullcase){
    FPR <- mean(pvals < 0.05)
    result <- list(FPR=FPR,
                   FWE=any(qvals < 0.05, na.rm=T))
  } else{

    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$exposure_eff != 0]
    DA_truth_binary <- taxa_truth$exposure_eff == 0
    roc_obj <-  suppressMessages(roc(DA_truth_binary, pvals))
    ancombc_DA_taxa <- output_df$taxon[output_df$diff_exposure1]
    check_DAtaxa <- ancombc_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


