suppressMessages(library(pROC))

# evaluate performance
evaluation <- function(gene_truth, ptda_result, nullcase=FALSE){
  ptda_pvals <- ptda_result$P_Value
  result <- NULL
  if (nullcase){
    result <- list(FPR = mean(ptda_pvals$pval < 0.05))
  } else{
    true_DA_gene <- rownames(gene_truth)[gene_truth$Logfold != 0]
    DA_truth_binary <- gene_truth[ptda_pvals$Gene, ] == 0
    roc_obj <- roc(DA_truth_binary, ptda_pvals$pval)
    check_refgene <- ptda_result$Reference_Gene %in% true_DA_gene
    check_DAgene <- ptda_result$DA_Gene %in% true_DA_gene
    refgene_error <- mean(check_refgene)
    FDR <- 1 - mean(check_DAgene)
    Power <- sum(check_DAgene) / length(true_DA_gene)
    result <- list(refgene_error=refgene_error,
                   FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  return(result)
}


