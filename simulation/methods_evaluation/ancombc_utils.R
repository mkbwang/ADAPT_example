suppressMessages(library(pROC))


# evaluate performance
evaluation_ancombc <- function(ancombc_result, taxa_truth=NULL, nullcase=FALSE){
  result <- NULL
  output_df <- ancombc_result$res
  pvals <- output_df$p_main1
  qvals <- output_df$q_main1
  if (nullcase){
    result <- list(FPR_1 = mean(pvals < 0.01),
                   FPR_5 = mean(pvals < 0.05),
                   FPR_10 = mean(pvals < 0.1),
                   num_FD = sum(qvals < 0.05))
  } else{

    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <-  suppressMessages(roc(DA_truth_binary, pvals))
    ancombc_DA_taxa <- output_df$taxon[output_df$diff_main1]
    check_DAtaxa <- ancombc_DA_taxa %in% true_DA_taxa
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


