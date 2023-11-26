suppressMessages(library(pROC))



# evaluate performance
evaluation_ancombc <- function(taxa_truth, ancombc_result, nullcase=FALSE){
  result <- NULL
  output_df <- ancombc_result$res
  pvals <- output_df$p_main1
  qvals <- output_df$q_main1
  if (nullcase){
    FPR <- mean(pvals < 0.05)
    result <- list(FPR=FPR,
                   FWE=any(pvals < 0.05/length(pvals), na.rm=T))
  } else{

    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <-  suppressMessages(roc(DA_truth_binary, pvals))
    ancombc_DA_taxa <- output_df$taxon[output_df$diff_main1]
    check_DAtaxa <- ancombc_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    FPR <- mean(pvals[DA_truth_binary] < 0.05)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   FPR=FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


