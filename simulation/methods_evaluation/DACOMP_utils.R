suppressMessages(library(pROC))

# evaluate performance
evaluation_dacomp <- function(dacomp_result, taxa_truth=NULL,  nullcase=FALSE){
  result <- NULL
  raw_pval <- dacomp_result$p.values.test
  adjusted_pval <- dacomp_result$p.values.test.adjusted
  if (nullcase){
    result <- list(
      FPR_1 = mean(raw_pval < 0.01, na.rm=TRUE),
      FPR_5 = mean(raw_pval < 0.05, na.rm=TRUE),
      FPR_10 = mean(raw_pval < 0.1, na.rm=TRUE),
      num_FD = sum(adjusted_pval < 0.05, na.rm=T))
  } else{
    
    dacomp_decision <- adjusted_pval < 0.05
    dacomp_decision[is.na(dacomp_decision)] <- FALSE
    
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <- suppressMessages(roc(DA_truth_binary[!is.na(raw_pval)], raw_pval[!is.na(raw_pval)]))

    check_DAtaxa <- dacomp_decision & !(DA_truth_binary)
    FDR <- 1 - sum(check_DAtaxa)/sum(dacomp_decision)
    FPR <- mean(raw_pval[DA_truth_binary] < 0.05, na.rm=T)
    Power <- sum(check_DAtaxa) / sum(!(DA_truth_binary))
    result <- list(FDR = FDR,
                   FPR = FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


