suppressMessages(library(pROC))

# evaluate performance
evaluation_dacomp <- function(taxa_truth, dacomp_result, nullcase=FALSE){
  result <- NULL
  raw_pval <- dacomp_result$p.values.test
  adjusted_pval <- dacomp_result$p.values.test.adjusted
  if (nullcase){
    FPR <- mean(raw_pval < 0.05, na.rm=TRUE)
    result <- list(FPR=FPR, FWE = any(adjusted_pval < 0.05, na.rm=T))
  } else{
    
    dacomp_decision <- adjusted_pval < 0.05
    dacomp_decision[is.na(dacomp_decision)] <- FALSE
    
    DA_truth_binary <- taxa_truth$exposure_eff == 0
    roc_obj <- suppressMessages(roc(DA_truth_binary[!is.na(raw_pval)], raw_pval[!is.na(raw_pval)]))

    check_DAtaxa <- dacomp_decision & !(DA_truth_binary)
    FDR <- 1 - sum(check_DAtaxa)/sum(dacomp_decision)
    Power <- sum(check_DAtaxa) / sum(!(DA_truth_binary))
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


