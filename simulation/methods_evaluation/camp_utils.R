
suppressMessages(library(pROC))



# evaluate performance
evaluation_camp <- function(camp_result, taxa_truth=NULL,  nullcase=FALSE){
  result <- NULL
  pvals <- camp_result[, 1]
  qvals <- p.adjust(pvals, method="BH")
  names(pvals) <- taxa_truth$taxaname
  if (nullcase){
    FPR <- mean(pvals < 0.05)
    result <- list(FPR_1 = mean(pvals < 0.01),
                   FPR_5 = mean(pvals < 0.05),
                   FPR_10 = mean(pvals < 0.1),
                   FWE_1 = any(pvals < 0.01/length(pvals), na.rm=T),
                   FWE_5 = any(pvals < 0.05/length(pvals), na.rm=T),
                   FWE_10 = any(pvals < 0.1/length(pvals), na.rm=T))
  } else{
    true_DA_taxa <- taxa_truth$taxaname[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <-  suppressMessages(roc(DA_truth_binary, pvals))
    camp_DA_taxa <- taxa_info$taxaname[qvals < 0.05]
    check_DAtaxa <- camp_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    FPR <- mean(pvals[DA_truth_binary] < 0.05)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   FPR=FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}




