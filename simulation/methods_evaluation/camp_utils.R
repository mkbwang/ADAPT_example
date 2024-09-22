
suppressMessages(library(pROC))



# evaluate performance
evaluation_camp <- function(camp_result, taxa_truth=NULL,  nullcase=FALSE){
  result <- NULL
  pvals <- camp_result[, 1]
  qvals <- p.adjust(pvals, method="BH")
  names(pvals) <- taxa_truth$taxaname
  if (nullcase){
    FPR <- mean(pvals < 0.05)
    result <- list(FPR_1 = mean(pvals < 0.01, na.rm=T),
                   FPR_5 = mean(pvals < 0.05, na.rm=T),
                   FPR_10 = mean(pvals < 0.1, na.rm=T),
                   num_FD = sum(qvals < 0.05, na.rm=T))
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




