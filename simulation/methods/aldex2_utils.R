suppressMessages(library(pROC))

# evaluate performance
evaluation_aldex2 <- function(taxa_truth, aldex_result, nullcase=FALSE){
  raw_pval <- aldex_result$`main:pval`
  adjusted_pval <- aldex_result$`main:pval.padj`
  result <- NULL
  if (nullcase){
    FPR <- mean(raw_pval < 0.05)
    result <- list(FPR=FPR, FWE=any(raw_pval < 0.05/nrow(aldex_result)))
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <- roc(DA_truth_binary, raw_pval)
    aldex_DA_taxa <- rownames(aldex_result)[adjusted_pval < 0.05]
    check_DAtaxa <- aldex_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    FPR <- mean(raw_pval[DA_truth_binary] < 0.05)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   FPR = FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }

  return(result)
}
