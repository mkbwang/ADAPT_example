
# overall settings for planned experiment

rm(list=ls())
folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation"

# Familywise error rates
FWER_settings <- expand.grid(Depth_confound=c(TRUE, FALSE), 
                             Nsample=c(50, 100, 200),
                             Zeroinflate=0,
                             Covar_confound=FALSE)

write.csv(FWER_settings, file.path(folder, 'FWER_settings.csv'), row.names=F)

# Proportion of DA taxa
PropDA_settings <- expand.grid(PropDA=c(0.05, 0.1, 0.2),
                               Direction=c("balanced", "unbalanced"),
                               FoldChange=4,
                               Nsample=100,
                               Depth=2e4,
                               Zeroinflate=0,
                               Depth_confound=FALSE,
                               Covar_confound=FALSE)


write.csv(PropDA_settings, file.path(folder, 'PropDA_settings.csv'), row.names=F)



# Fold Changes
FoldChange_settings <- expand.grid(PropDA=0.1,
                                   Direction=c("balanced", "unbalanced"),
                                   FoldChange=c(3,4,5),
                                   Nsample=100,
                                   Depth=2e4,
                                   Zeroinflate=0,
                                   Depth_confound=FALSE,
                                   Covar_confound=FALSE)

write.csv(FoldChange_settings, file.path(folder, 'FoldChange_settings.csv'), row.names=F)

# Sample Sizes
SampleSize_settings <- expand.grid(PropDA=0.1,
                                   Direction=c("balanced", "unbalanced"),
                                   FoldChange=4,
                                   Nsample=c(50, 100, 200),
                                   Depth=2e4,
                                   Zeroinflate=0,
                                   Depth_confound=FALSE,
                                   Covar_confound=FALSE)

write.csv(SampleSize_settings, file.path(folder, 'SampleSize_settings.csv'), row.names=F)



# Sequence Depths
SeqDepth_settings <- expand.grid(PropDA=0.1,
                                 Direction=c("balanced", "unbalanced"),
                                 FoldChange=4,
                                 Nsample=100,
                                 Depth=c(1e4, 2e4, 5e4),
                                 Zeroinflate=0,
                                 Depth_confound=FALSE,
                                 Covar_confound=FALSE)


write.csv(SeqDepth_settings, file.path(folder, 'SeqDepth_settings.csv'), row.names=F)



# Zero Inflation
ZeroInflation_settings <- expand.grid(PropDA=0.1,
                                      Direction=c("balanced", "unbalanced"),
                                      FoldChange=4,
                                      Nsample=100,
                                      Depth=2e4,
                                      Zeroinflate=c(0, 0.2, 0.4),
                                      Depth_confound=FALSE,
                                      Covar_confound=FALSE)


write.csv(ZeroInflation_settings, file.path(folder, 'ZeroInflation_settings.csv'), row.names=F)



