suppressMessages(library(pROC))


# evaluate performance
evaluation_linda <- function(linda_result, taxa_truth=NULL, nullcase=FALSE){
  result <- NULL
  linda_df <- linda_result$output$main
  raw_pval <- linda_df$pvalue
  adjusted_pval <- linda_df$padj
  if (nullcase){
    FPR <- mean(raw_pval < 0.05, na.rm=TRUE)
    result <- list(
      FPR_1 = mean(raw_pval < 0.01, na.rm=TRUE),
      FPR_5 = mean(raw_pval < 0.05, na.rm=TRUE),
      FPR_10 = mean(raw_pval < 0.1, na.rm=TRUE),
      num_FD = sum(adjusted_pval < 0.05))
  } else{
    linda_decision <- linda_df$reject
    
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <- suppressMessages(roc(DA_truth_binary[!is.na(raw_pval)], raw_pval[!is.na(raw_pval)]))

    check_DAtaxa <- linda_decision & !(DA_truth_binary)
    FDR <- 1 - sum(check_DAtaxa)/sum(linda_decision)
    FPR <- mean(raw_pval[DA_truth_binary] < 0.05)
    Power <- sum(check_DAtaxa) / sum(!(DA_truth_binary))
    result <- list(FDR = FDR,
                   FPR = FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


