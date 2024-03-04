suppressMessages(library(pROC))

# evaluate performance
evaluation_zicoseq <- function(zicoseq_result, taxa_truth=NULL, nullcase=FALSE){
  result <- NULL
  raw_pval <- zicoseq_result$p.raw
  adjusted_pval <- zicoseq_result$p.adj.fdr
  fwer_pval <- zicoseq_result$p.adj.fwer
  if (nullcase){
    FPR <- mean(raw_pval < 0.05, na.rm=TRUE)
    result <- list(FPR_1 = mean(raw_pval < 0.01, na.rm=TRUE), 
                   FPR_5 = mean(raw_pval < 0.05, na.rm=TRUE),
                   FPR_10 = mean(raw_pval < 0.1, na.rm=TRUE), 
                   FWE_1 = any(fwer_pval < 0.01, na.rm=T),
                   FWE_5 = any(fwer_pval < 0.05, na.rm=T),
                   FWE_10 = any(fwer_pval < 0.1, na.rm=T))
  } else{
    
    zicoseq_decision <- adjusted_pval < 0.05
    zicoseq_decision[is.na(zicoseq_decision)] <- FALSE
    
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <- suppressMessages(roc(DA_truth_binary[!is.na(raw_pval)], raw_pval[!is.na(raw_pval)]))

    check_DAtaxa <- zicoseq_decision & !(DA_truth_binary)
    FDR <- 1 - sum(check_DAtaxa)/sum(zicoseq_decision)
    FPR <- mean(raw_pval[DA_truth_binary] < 0.05)
    Power <- sum(check_DAtaxa) / sum(!(DA_truth_binary))
    result <- list(FDR = FDR,
                   FPR = FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


