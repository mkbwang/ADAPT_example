
rm(list=ls())

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation"
suppressMessages(source(file.path(folder, "SparseDOSSA_setting.R")))
# suppressMessages(source(file.path(folder, "camp.R")))

baseline_params <- readRDS(file.path(folder, 'sparsedossa_baseline.rds'))

settings_df <- expand.grid(nSample=100,
                           nTaxa=500,
                           seqdepth_mean=2e4,
                           propDA=c(0.05, 0.1, 0.2, 0.3), # c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
                           main_mean_logfold=log(5),
                           direction=c("balanced", "unbalanced"),
                           propcf=0,
                           cf_sd=1,
                           cf_main_corr=0,
                           cf_mean_logfold=1,
                           stringsAsFactors = FALSE)

# parse the arguments
library(optparse)
option_list <- list(make_option(c("-c", "--choice"), type="integer", default=1, 
                                help="Choice of Setting [default=1]"),
                    make_option(c("-s", "--seed"), type="integer", default=1, 
                                help="seed [default=1]"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
choice <- opt$choice
myseed <- opt$seed


sparsedossa_setting <- SimulateSettings(baseline=baseline_params,
                                        nSample=settings_df$nSample[choice],
                                        nTaxa=settings_df$nTaxa[choice],
                                        seqdepth_logmean=log(settings_df$seqdepth_mean[choice]),
                                        propDA=settings_df$propDA[choice],
                                        main_mean_logfold=settings_df$main_mean_logfold[choice],
                                        direction=settings_df$direction[choice],
                                        propcf=settings_df$propcf[choice],
                                        cf_sd=settings_df$cf_sd[choice],
                                        cf_main_corr=settings_df$cf_main_corr[choice],
                                        cf_mean_logfold = settings_df$cf_mean_logfold[choice],
                                        seed=myseed)

metadata_df <- data.frame(sparsedossa_setting$sample_metadata)
colnames(metadata_df) <- c("main", "confounder")
rownames(metadata_df) <- sprintf("Sample%d", seq(1, nrow(metadata_df)))


simulated_data <- SparseDOSSA2(template=sparsedossa_setting$baseline_param,
                               new_features=TRUE,
                               n_sample = nrow(sparsedossa_setting$sample_metadata),
                               n_feature=nrow(sparsedossa_setting$taxa_info),
                               spike_metadata=sparsedossa_setting$spike_metadata,
                               metadata_matrix=sparsedossa_setting$sample_metadata,
                               median_read_depth=exp(sparsedossa_setting$baseline_param$depth_fit[1]),
                               verbose=F)


count_table <- simulated_data$simulated_data
rownames(count_table) <- sprintf("Taxon%d", seq(1, nrow(sparsedossa_setting$taxa_info)))
taxa_info <- sparsedossa_setting$taxa_info


prevalence_cutoff <- 0.05
seqdepth_cutoff <- 1e3
prevalences <- rowMeans(count_table > 0)
taxa_filter <- prevalences > 0.05
count_table <- count_table[taxa_filter, ]
taxa_info <- taxa_info[taxa_filter, , drop=F]
true_DAtaxa <- taxa_info$taxaname[taxa_info$log_main_eff != 0]

seqdepths <- colSums(count_table)


seqdepth_filter <- seqdepths > seqdepth_cutoff
count_table <- count_table[, seqdepth_filter, drop=F]
metadata_df <- metadata_df[seqdepth_filter, ]



suppressMessages(library(phyloseq))
phyobj <- phyloseq(otu_table(count_table, taxa_are_rows = T),
                   sample_data(metadata_df))


# ADAPT
suppressMessages(library(ADAPT))
begin <- proc.time()
adapt_output <- adapt(input_data = phyobj, cond.var="main")
adapt_time <- proc.time() - begin
adapt_duration_noboot <- adapt_time[3]
source(file.path(folder, "methods_evaluation", "adapt_utils.R"))
adapt_performance <- suppressMessages(evaluation_adapt(taxa_truth=taxa_info,
                                                       adapt_result=adapt_output,
                                                       nullcase=F))




# use all the nonDA taxa as reference

nonDA_taxa <- taxa_info$taxaname[taxa_info$log_main_eff == 0]
alternative_result <- ADAPT:::count_ratio(count_table = t(count_table), design_matrix=cbind(1, metadata_df$main),
                                          reftaxa = nonDA_taxa, censor=1, test_all=TRUE)
alternative_result$adjpval <- p.adjust(alternative_result$pval, method="BH")
alternative_DAtaxa <- alternative_result$Taxa[alternative_result$adjpval < 0.05]


alternative_FDR <- 1-mean(alternative_DAtaxa %in% true_DAtaxa)
alternative_power <- mean(true_DAtaxa %in% alternative_DAtaxa)
num_intersection <- length(intersect(alternative_DAtaxa, adapt_output@signal))
num_union <- length(union(alternative_DAtaxa, adapt_output@signal))
IoU <- length(num_intersection) / length(num_union)

detailed_performance <- list(num_taxa = adapt_performance$num_taxa, num_reftaxa=adapt_performance$num_reftaxa, 
                             ref_error = adapt_performance$ref_error,  num_DAtaxa=length(adapt_output@signal), 
                             FDR = adapt_performance$FDR, 
                             power = adapt_performance$Power,
                             num_altDAtaxa = length(alternative_DAtaxa),
                             altFDR = alternative_FDR,
                             altPower = alternative_power,
                             num_intersection = num_intersection,
                             num_union = num_union,
                             IoU = IoU)


output_file <- sprintf("SparseDOSSA_%d_%d_adapt_details.rds", choice, myseed)
saveRDS(detailed_performance, file=file.path(folder, "propDA", "details_ADAPT", output_file))


