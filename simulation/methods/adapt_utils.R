suppressMessages(library(pROC))

# evaluate performance
evaluation_adapt <- function(adapt_result, taxa_truth=NULL, nullcase=FALSE){
  
  adapt_pvals <- adapt_result$P_Value
  adapt_DAtaxa <- adapt_result$DA_Taxa
  adapt_reference <- adapt_result$Reference_Taxa
  result <- NULL
  if (nullcase){
    result <- list(
      FPR_1 = mean(adapt_pvals$pval < 0.01, na.rm=T),
      FPR_5 = mean(adapt_pvals$pval < 0.05, na.rm=T),
      FPR_10 = mean(adapt_pvals$pval < 0.1, na.rm=T),
      FWE_1 = any(adapt_pvals$pval < 0.01 / length(adapt_pvals$pval), na.rm=T),
      FWE_5 = any(adapt_pvals$pval < 0.05 / length(adapt_pvals$pval), na.rm=T),
      FWE_10 = any(adapt_pvals$pval < 0.1 / length(adapt_pvals$pval), na.rm=T))
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <- roc(DA_truth_binary, adapt_pvals$pval)
    check_ref <- adapt_reference %in% true_DA_taxa
    check_DA <- adapt_DAtaxa %in% true_DA_taxa
    ref_error <- mean(check_ref)
    FDR <- 1 - mean(check_DA)
    FPR <- mean(adapt_pvals$pval[DA_truth_binary] < 0.05)
    Power <- sum(check_DA) / length(true_DA_taxa)
    result <- list(ref_error=ref_error,
                   FDR = FDR,
                   FPR = FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  return(result)
}


