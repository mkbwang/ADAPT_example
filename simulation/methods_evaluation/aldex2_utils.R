suppressMessages(library(pROC))

# evaluate performance
evaluation_aldex2<- function( aldex_result, test=c("ttest", "glm"), taxa_truth=NULL, nullcase=FALSE){
  if (test == "ttest"){
    raw_pval <- aldex_result$wi.ep
    adjusted_pval <- aldex_result$wi.eBH
  } else{
    raw_pval <- aldex_result$`main:pval`
    adjusted_pval <- aldex_result$`main:pval.padj`
  }
  result <- NULL
  if (nullcase){
    result <- list(
      FPR_1 = mean(raw_pval < 0.01),
      FPR_5 = mean(raw_pval < 0.05),
      FPR_10 = mean(raw_pval < 0.1),
      FWE_1 =any(raw_pval < 0.01/nrow(aldex_result)),
      FWE_5 =any(raw_pval < 0.05/nrow(aldex_result)),
      FWE_10 =any(raw_pval < 0.1/nrow(aldex_result)))
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

