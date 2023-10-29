
rm(list=ls())

# Familywise error rates
FWER_settings <- expand.grid(Depth_confound=c(TRUE, FALSE), 
                             Nsample=c(50, 100, 200),
                             Zeroinflate=0,
                             Depth_confound=FALSE,
                             Covar_confound=FALSE)

# Proportion of DA taxa
PropDA_settings <- expand.grid(PropDA=c(0.05, 0.1, 0.2),
                               Direction=c("balanced", "unbalanced"),
                               FoldChange=4,
                               Nsample=100,
                               Depth=2e4,
                               Zeroinflate=0,
                               Depth_confound=FALSE,
                               Covar_confound=FALSE)


# Fold Changes
FoldChange_settings <- expand.grid(PropDA=0.1,
                                   Direction=c("balanced", "unbalanced"),
                                   FoldChange=c(3,4,5),
                                   Nsample=100,
                                   Depth=2e4,
                                   Zeroinflate=0,
                                   Depth_confound=FALSE,
                                   Covar_confound=FALSE)


# Sample Sizes
SampleSize_settings <- expand.grid(PropDA=0.1,
                                   Direction=c("balanced", "unbalanced"),
                                   FoldChange=4,
                                   Nsample=c(50, 100, 200),
                                   Depth=2e4,
                                   Zeroinflate=0,
                                   Depth_confound=FALSE,
                                   Covar_confound=FALSE)



# Sequence Depths
SeqDepth_settings <- expand.grid(PropDA=0.1,
                                 Direction=c("balanced", "unbalanced"),
                                 FoldChange=4,
                                 Nsample=100,
                                 Depth=c(1e4, 2e4, 5e4),
                                 Zeroinflate=0,
                                 Depth_confound=FALSE,
                                 Covar_confound=FALSE)


# Zero Inflation
ZeroInflation_settings <- expand.grid(PropDA=0.1,
                                      Direction=c("balanced", "unbalanced"),
                                      FoldChange=4,
                                      Nsample=100,
                                      Depth=c(1e4, 2e4, 5e4),
                                      Zeroinflate=0,
                                      Depth_confound=FALSE,
                                      Covar_confound=FALSE)



# Depth Confounding
DepthConfound_settings <- expand.grid(PropDA=0.1,
                                      Direction=c("balanced", "unbalanced"),
                                      FoldChange=4,
                                      Nsample=100,
                                      Depth=2e4,
                                      Zeroinflate=0,
                                      Depth_confound=c(FALSE, TRUE),
                                      Covar_confound=FALSE)


