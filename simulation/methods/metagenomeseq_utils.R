suppressMessages(library(pROC))


# evaluate performance
evaluation_metagenomeseq <- function(taxa_truth, MR_result, nullcase=FALSE){
  
  performance <- NULL
  MR_coefficients <- MRcoefs(MR_result, number=nrow(taxa_truth), adjustMethod="BH")
  
  result <- NULL
  if (nullcase){
    FPR <- mean(MR_coefficients$pvalues < 0.05,  na.rm=TRUE)
    result <- list(FPR=FPR,
                   FWE=any(MR_coefficients$adjPvalues < 0.05, na.rm=T))
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$exposure_eff != 0]
    DA_truth_binary <- taxa_truth[rownames(MR_coefficients), 'exposure_eff'] == 0
    roc_obj <- roc(DA_truth_binary, MR_coefficients$pvalues, na.rm=TRUE)
    MR_DA_taxa <- rownames(MR_coefficients)[MR_coefficients$adjPvalues < 0.05]
    check_DAtaxa <- MR_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}
