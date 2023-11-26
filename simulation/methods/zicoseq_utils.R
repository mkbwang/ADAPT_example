suppressMessages(library(pROC))

# evaluate performance
evaluation_zicoseq <- function(taxa_truth, zicoseq_result, nullcase=FALSE){
  result <- NULL
  raw_pval <- zicoseq_result$p.raw
  adjusted_pval <- zicoseq_result$p.adj.fdr
  fwer_pval <- zicoseq_result$p.adj.fwer
  if (nullcase){
    FPR <- mean(raw_pval < 0.05, na.rm=TRUE)
    result <- list(FPR=FPR, FWE = any(fwer_pval < 0.05, na.rm=T))
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


