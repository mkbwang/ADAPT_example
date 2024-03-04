suppressMessages(library(pROC))



# evaluate performance
evaluation_maaslin2 <- function(maaslin2_result, taxa_truth=NULL,  nullcase=FALSE){
  result <- NULL
  maaslin2_result <- maaslin2_result[maaslin2_result$metadata == "main", ]
  pvals <- maaslin2_result$pval
  names(pvals) <- maaslin2_result$feature
  if (!nullcase){
    pvals <- pvals[rownames(taxa_truth)] 
  }
  qvals <- maaslin2_result$qval
  if (!nullcase){
    names(qvals) <- maaslin2_result$feature 
  }
  if (nullcase){
    FPR <- mean(pvals < 0.05)
    result <- list(FPR_1 = mean(pvals < 0.01),
                   FPR_5 = mean(pvals < 0.05),
                   FPR_10 = mean(pvals < 0.1),
                   FWE_1 = any(pvals < 0.01/length(pvals), na.rm=T),
                   FWE_5 = any(pvals < 0.05/length(pvals), na.rm=T),
                   FWE_10 = any(pvals < 0.1/length(pvals), na.rm=T))
  } else{
    
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$log_main_eff != 0]
    DA_truth_binary <- taxa_truth$log_main_eff == 0
    roc_obj <-  suppressMessages(roc(DA_truth_binary, pvals))
    maaslin2_DA_taxa <- maaslin2_result$feature[maaslin2_result$qval < 0.05]
    check_DAtaxa <- maaslin2_DA_taxa %in% true_DA_taxa
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


