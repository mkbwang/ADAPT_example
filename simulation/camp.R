#' Censoring-based analysis of microbiome proportions (CAMP)
#'
#' X A size (n x p) matrix, p is the number of taxa
#' cov A length n vector (in case of only group label) or a size (n x q) data frame, q is the number of covariates (including group)
#' returns A length p vector of p-values, or a size (p x q) data frame of p-values corresponding the covariates
library(survival)
camp <- function(X, cov) {
  n <- dim(X)[1] # Sample size
  p <- dim(X)[2] # Number of taxa
  Delta <- (X != 0) + 0 # Creating censoring indicator for survival analysis
  T <- X
  T[which(T == 0, arr.ind = TRUE)] <- min(X[X != 0]) # Zero replacement
  
  # Normalize the data
  T <- T / matrix(rowSums(T), n, p)
  T <- -log(T) # "Time" data
  
  if (is.vector(cov)) {
    pval.vec <- rep(NA, p)
    for (i in 1:p) {
      fit <- survdiff(Surv(T[,i], Delta[,i]) ~ cov)
      pval.vec[i] <- 1 - pchisq(fit$chisq, 1)
    }
    return(pval.vec)
  } else if (is.data.frame(cov)) {
    q <- dim(cov)[2]
    pval.df <- data.frame(matrix(NA, p, q))
    colnames(pval.df) <- names(cov)
    for (i in 1:p) {
      formula_str <- paste0("Surv(T[,", i, "], Delta[,", i, "]) ~ ", paste(names(cov), collapse = " + "))
      fit <- coxph(as.formula(formula_str), data=cov)
      pval.df[i, ] <- summary(fit)$coefficients[1:q, 5]
    }
    return(pval.df)
  } else {
    print("Cov is of the wrong data type!")
  }
}