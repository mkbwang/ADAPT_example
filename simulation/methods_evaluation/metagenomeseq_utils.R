suppressMessages(library(pROC))


# evaluate performance
evaluation_metagenomeseq <- function(MR_result, taxa_truth=NULL,  nullcase=FALSE){
  
  performance <- NULL
  MR_coefficients <- MRcoefs(MR_result, number=nrow(taxa_truth), adjustMethod="BH")
  
  result <- NULL
  if (nullcase){
    result <- list(FPR_1 = mean(MR_coefficients$pvalues < 0.01,  na.rm=TRUE),
                   FPR_5 = mean(MR_coefficients$pvalues < 0.05,  na.rm=TRUE),
                   FPR_10 = mean(MR_coefficients$pvalues < 0.1,  na.rm=TRUE),
                   num_FD = sum(MR_coefficients$adjPvalues < 0.05, na.rm=T))
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth[rownames(MR_coefficients), 'log_main_eff'] == 0
    roc_obj <- roc(DA_truth_binary, MR_coefficients$pvalues, na.rm=TRUE)
    MR_DA_taxa <- rownames(MR_coefficients)[MR_coefficients$adjPvalues < 0.05]
    check_DAtaxa <- MR_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    FPR <- mean(MR_coefficients$pvalues[DA_truth_binary] < 0.05, na.rm=T)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   FPR = FPR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}
