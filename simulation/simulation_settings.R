
# overall settings for planned experiment

rm(list=ls())
folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation"

# T1errors
settings_df <- expand.grid(nSample=c("Sample Size = 50", "Sample Size = 100"),
                           nTaxa=500,
                           depth_fold=c("Balanced Library Size", "Unbalanced Library Size (4-fold)",
                                        "Unbalanced Library Size (10-fold)"))
write.csv(settings_df, file.path(folder, 'T1error_settings.csv'), row.names=F)


# Fold Changes
settings_df <- expand.grid(FoldChange=c(3, 3.5, 4, 5, 6),
                           direction=c("Balanced Change", "Unbalanced Change"),
                           stringsAsFactors = FALSE)
write.csv(settings_df, file.path(folder, 'FoldChange_settings.csv'), row.names=F)



# Proportion of DA taxa
settings_df <- expand.grid(propDA=c(0.05, 0.1, 0.2, 0.3),
                           direction=c("Balanced Change", "Unbalanced Change"),
                           stringsAsFactors = FALSE)
write.csv(settings_df, file.path(folder, 'PropDA_settings.csv'), row.names=F)



# Sample Sizes
settings_df <- expand.grid(nSample=c(50, 80, 100, 150, 200),
                           direction=c("Balanced Change", "Unbalanced Change"),
                           stringsAsFactors = FALSE)
write.csv(settings_df, file.path(folder, 'SampleSize_settings.csv'), row.names=F)



# Sequence Depths
settings_df <- expand.grid(seqdepth_mean=c(5e3, 1e4, 2e4, 4e4),
                           direction=c("Balanced Change", "Unbalanced Change"),
                           stringsAsFactors = FALSE)
write.csv(settings_df, file.path(folder, 'SeqDepth_settings.csv'), row.names=F)


# Taxa Number
settings_df <- expand.grid(nTaxa=c(100, 200, 500, 1000),
                           direction=c("Balanced Change", "Unbalanced Change"),
                           stringsAsFactors = FALSE)
write.csv(settings_df, file.path(folder, 'TaxaNum_settings.csv'), row.names=F)


# confounder
settings_df <- expand.grid(propcf=c(0, 0.2, 0.4),
                           cf_main_corr=c(0, 0.3, 0.6),
                           stringsAsFactors = FALSE)
write.csv(settings_df, file.path(folder, 'confounder_settings.csv'), row.names=F)



