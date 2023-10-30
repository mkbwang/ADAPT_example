suppressMessages(library(pROC))


# evaluate performance
evaluation_linda <- function(taxa_truth, linda_result, nullcase=FALSE){
  result <- NULL
  linda_df <- linda_result$output$exposure
  raw_pval <- linda_df$pvalue
  adjusted_pval <- linda_df$padj
  if (nullcase){
    FPR <- mean(raw_pval < 0.05, na.rm=TRUE)
    result <- list(FPR=FPR, FWE = any(linda_df$reject, na.rm=T))
  } else{
    linda_decision <- linda_df$reject
    
    DA_truth_binary <- taxa_truth$exposure_eff == 0
    roc_obj <- suppressMessages(roc(DA_truth_binary[!is.na(raw_pval)], raw_pval[!is.na(raw_pval)]))

    check_DAtaxa <- linda_decision & !(DA_truth_binary)
    FDR <- 1 - sum(check_DAtaxa)/sum(linda_decision)
    Power <- sum(check_DAtaxa) / sum(!(DA_truth_binary))
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


