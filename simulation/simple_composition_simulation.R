

library(MASS)

library(extraDistr)

group1_abs <- rpois() 
foldchange <- exp(rnorm(100, mean=0, sd=1))

alldata <- matrix()

group1 <- matrix(rnegbin(20000, mu=1000, theta=5), nrow=100) 
group2_same <- matrix(rnegbin(10000, mu=1000, theta=5), nrow=100)
group2_decrease <- matrix(rnegbin(5000, mu=500, theta=5), nrow=100)
group2_increase <- matrix(rnegbin(5000, mu=4000, theta=5), nrow=100)
group2 <- cbind(group2_same, group2_decrease, group2_increase)

alldata <- rbind(group1, group2)
colnames(alldata) <- sprintf("Taxon_%d", seq(1, 200))


sequencing_depths <- rnegbin(n=200, mu=10000, theta=5)
count_mat <- matrix(0,  nrow=200, ncol=200)
colnames(count_mat) <- sprintf("Taxon_%d", seq(1, 200))

for (j in 1:nrow(count_mat)){
  count_mat[j,] <- rmvhyper(nn=1, n=alldata[j, ], k=sequencing_depths[j])
}

logit_result <- data.frame(Taxa_names = colnames(alldata),
                           Logit_effect_raw = rep(0, nrow(alldata)),
                           Logit_SE_raw = rep(0, nrow(alldata)),
                           Logit_effect_corrected = rep(0, nrow(alldata)),
                           Logit_SE_corrected = rep(0, nrow(alldata)),
                           Logit_Pval = rep(0, nrow(alldata)))

covariate <- c(rep(0, 100), rep(1, 100))

for (j in 1:ncol(count_mat)){
  taxon_interest <- count_mat[, j]
  taxa_others <- rowSums(count_mat[, -j])
  test_df <- cbind(taxon_interest, taxa_others, covariate) |> as.data.frame()
  colnames(test_df) <- c("Taxon", "Others", 'CRC')
  rawglm <- glm(cbind(Taxon, Others) ~ CRC, family=binomial(link = "logit"),
                data=test_df)
  rawglm_summary <- summary(rawglm)
  corrected_glm <- glm.binomial.disp(rawglm, verbose = FALSE) %>%
    summary()
  logit_result$Logit_effect_raw[j] <- rawglm_summary$coefficients[2, 1]
  logit_result$Logit_SE_raw[j] <- rawglm_summary$coefficients[2, 2]
  logit_result$Logit_effect_corrected[j] <- corrected_glm$coefficients[2, 1]
  logit_result$Logit_SE_corrected[j] <- corrected_glm$coefficients[2, 2]
  logit_result$Logit_Pval[j] <- corrected_glm$coefficients[2, 4]
}

