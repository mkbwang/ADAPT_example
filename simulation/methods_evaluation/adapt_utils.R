suppressMessages(library(pROC))

# evaluate performance
evaluation_adapt <- function(adapt_result, taxa_truth=NULL, nullcase=FALSE){
  
  adapt_pvals <- adapt_result@details$pval
  adapt_DAtaxa <- adapt_result@signal
  adapt_reference <- adapt_result@reference
  result <- NULL
  if (nullcase){
    result <- list(
      FPR_1 = mean(adapt_pvals < 0.01, na.rm=T),
      FPR_5 = mean(adapt_pvals < 0.05, na.rm=T),
      FPR_10 = mean(adapt_pvals < 0.1, na.rm=T),
      FWE_1 = any(adapt_pvals < 0.01 / length(adapt_pvals), na.rm=T),
      FWE_5 = any(adapt_pvals < 0.05 / length(adapt_pvals), na.rm=T),
      FWE_10 = any(adapt_pvals < 0.1 / length(adapt_pvals), na.rm=T))
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <- roc(DA_truth_binary, adapt_pvals)
    check_ref <- adapt_reference %in% true_DA_taxa
    check_DA <- adapt_DAtaxa %in% true_DA_taxa
    ref_error <- sum(check_ref)
    FDR <- 1 - mean(check_DA)
    FPR <- mean(adapt_pvals[DA_truth_binary] < 0.05)
    Power <- sum(check_DA) / length(true_DA_taxa)
    result <- list(num_taxa = length(adapt_pvals),
                   num_DAtaxa = sum(!DA_truth_binary),
                   num_reftaxa = length(adapt_reference),
                   ref_error=ref_error,
                   FDR = FDR,
                   FPR = FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  return(result)
}


