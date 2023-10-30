suppressMessages(library(pROC))

# evaluate performance
evaluation_adapt <- function(taxa_truth, adapt_result, nullcase=FALSE){
  
  adapt_pvals <- adapt_result$P_Value
  adapt_DAtaxa <- adapt_result$DA_Gene
  adapt_reference <- adapt_result$Reference_Gene
  result <- NULL
  if (nullcase){
    result <- list(FPR = mean(adapt_pvals$pval < 0.05),
                   FWE = length(adapt_DAtaxa) > 0)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$exposure_eff != 0]
    DA_truth_binary <- taxa_truth$exposure_eff == 0
    roc_obj <- roc(DA_truth_binary, adapt_pvals$pval)
    check_ref <- adapt_reference %in% true_DA_taxa
    check_DA <- adapt_DAtaxa %in% true_DA_taxa
    ref_error <- mean(check_ref)
    FDR <- 1 - mean(check_DA)
    Power <- sum(check_DA) / length(true_DA_taxa)
    result <- list(ref_error=ref_error,
                   FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  return(result)
}


