
rm(list=ls())

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation"
suppressMessages(source(file.path(folder, "SparseDOSSA_setting.R")))
suppressMessages(source(file.path(folder, "camp.R")))

baseline_params <- readRDS(file.path(folder, 'sparsedossa_baseline.rds'))

settings_df <- expand.grid(nSample=100,
                           nTaxa=500,
                           seqdepth_mean=2e4,
                           propDA=0.1,
                           main_mean_logfold=c(log(3), log(3.5), log(4), log(5), log(6)),
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
adapt_output <- adapt(input_data=phyobj, cond.var="main")
adapt_time <- proc.time() - begin
adapt_duration_noboot <- adapt_time[3]
source(file.path(folder, "methods_evaluation", "adapt_utils.R"))
adapt_performance <- suppressMessages(evaluation_adapt(taxa_truth=taxa_info,
                                                       adapt_result=adapt_output,
                                                       nullcase=F))



# Aldex2
suppressMessages(library(ALDEx2))
begin <- proc.time()
aldex_result <- suppressMessages(aldex(reads=count_table, conditions=metadata_df$main, mc.samples=128,
                                       test="t"))
aldex_time <- proc.time() - begin
aldex_duration <- aldex_time[3]
source(file.path(folder, "methods_evaluation", "aldex2_utils.R"))
aldex_performance <- suppressMessages(evaluation_aldex2(taxa_truth=taxa_info, 
                                       aldex_result=aldex_result, test="ttest",
                                       nullcase = F))


# Maaslin2
suppressMessages(library(Maaslin2))
count_df <- data.frame(t(count_table))
begin <- proc.time()
maaslin2_output <-invisible(Maaslin2(input_data=count_df, input_metadata=metadata_df,
                            output = file.path(folder, 'FoldChange', 'maaslin2_output', sprintf("seed_%d", myseed)),
                            min_prevalence=0.05,
                            max_significance = 0.05,
                            fixed_effects=c("main", "confounder"),
                            cores=4, plot_heatmap = F,plot_scatter = F))
maaslin2_time <- proc.time() - begin
maaslin2_duration <- maaslin2_time[3]
maaslin2_result <- maaslin2_output$results
source(file.path(folder, "methods_evaluation", "maaslin2_utils.R"))
maaslin2_performance <- suppressMessages(evaluation_maaslin2(taxa_truth=taxa_info,
                                                             maaslin2_result=maaslin2_result,
                                                             nullcase=F))


# MetagenomeSeq
suppressMessages(library(metagenomeSeq))
MRmetadata <- AnnotatedDataFrame(metadata_df)
MRobj <- newMRexperiment(counts=count_table,
                         phenoData=MRmetadata)
begin <- proc.time()
MRobj <- wrenchNorm(MRobj, condition=metadata_df$main)
mod <- model.matrix(~main, data=pData(MRobj))
MRfit <- fitFeatureModel(MRobj, mod)
metagenomeseq_time <- proc.time() - begin
metagenomeseq_duration <- metagenomeseq_time[3]
source(file.path(folder, "methods_evaluation", "metagenomeseq_utils.R"))
metagenomeseq_performance <- suppressMessages(evaluation_metagenomeseq(taxa_truth=taxa_info,
                                           MR_result=MRfit, nullcase=F))



# DACOMP
library(dacomp)
begin <- proc.time()
selected_references <- dacomp.select_references(
  X = t(count_table), verbose = F, run_in_parallel = T, Nr.Cores=4,
  minimal_TA=10)
updated_selected_references <- dacomp.validate_references(X = t(count_table), #Counts
                                                          Y = metadata_df$main, #Traits
                                                          ref_obj = selected_references, #reference checked, must be object from dacomp.select_references(...)
                                                          test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                          Q_validation = 0.1, #FDR level for checking
                                                          Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                          Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                          NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                          Verbose = F)
dacomp_output <- dacomp.test(X = t(count_table), #counts data
                          y = metadata_df$main, #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = updated_selected_references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          disable_DSFDR = T,
                          q=0.05,
                          verbose = F)
dacomp_time <- proc.time() - begin
dacomp_duration <- dacomp_time[3]
source(file.path(folder, "methods_evaluation", "DACOMP_utils.R"))
dacomp_performance <- evaluation_dacomp(taxa_truth=taxa_info,
                                        dacomp_result=dacomp_output,
                                        nullcase=F)



# ZicoSeq
library(GUniFrac)
begin <- proc.time()
zicoseq_output <- ZicoSeq(meta.dat=metadata_df,
                          feature.dat=count_table,
                          grp.name="main",
                          adj.name="confounder",is.fwer=T,
                          perm.no=999)
zicoseq_result <- cbind.data.frame(p.raw = zicoseq_output$p.raw,
                                p.adj.fdr = zicoseq_output$p.adj.fdr,
                                p.adj.fwer = zicoseq_output$p.adj.fwer)
zicoseq_time <- proc.time()-begin
zicoseq_duration <- zicoseq_time[3]
source(file.path(folder, "methods_evaluation", "zicoseq_utils.R"))
zicoseq_performance <- evaluation_zicoseq(taxa_truth=taxa_info,
                                          zicoseq_result=zicoseq_result,
                                          nullcase=F)


#ANCOM
suppressMessages(library(ANCOMBC))

begin <- proc.time()
ancom_output <- ancom(data=phyobj, p_adj_method="BH", prv_cut=0.05,
                      main_var="main", adj_formula="confounder", n_cl=4)
ancom_time <- proc.time()-begin
ancom_duration <- ancom_time[3]
source(file.path(folder, "methods_evaluation", "ancom_utils.R"))
ancom_performance <- evaluation_ancom(taxa_truth=taxa_info,
                                      ancom_result=ancom_output$res,
                                      nullcase=F)


# ANCOMBC
ptm <- proc.time()
ancombc_output <- ancombc2(phyobj, fix_formula = 'main + confounder', p_adj_method='BH', global=T,
                          group='main', struc_zero=FALSE, prv_cut=0.05, n_cl=4)
ancombc_time <- proc.time() - ptm
ancombc_duration <- ancombc_time[3]
source(file.path(folder, "methods_evaluation", "ancombc_utils.R"))
ancombc_performance <- evaluation_ancombc(taxa_truth=taxa_info,
                                          ancombc_result=ancombc_output,
                                          nullcase=F)


# LinDA
suppressMessages(library(MicrobiomeStat))
begin <- proc.time()
linda_output <- MicrobiomeStat::linda(feature.dat=count_table,
                      meta.dat=metadata_df,
                      formula="~main+confounder",
                      n.cores=4)
linda_time <- proc.time() - begin
linda_duration <- linda_time[3]
source(file.path(folder, "methods_evaluation", "linda_utils.R"))
linda_performance <- evaluation_linda(taxa_truth=taxa_info,
                                      linda_result=linda_output,
                                      nullcase=F)

# CAMP
begin <- proc.time()
camp_output <- camp(X=t(count_table), cov=metadata_df)
camp_time <- proc.time() - begin
camp_duration <- camp_time[3]
source(file.path(folder, "methods_evaluation", "camp_utils.R"))
camp_performance <- evaluation_camp(camp_result = camp_output,
                                    taxa_truth = taxa_info,
                                    nullcase=FALSE)


# aggregate all the performances
Methods <- c("ADAPT", "ALDEx2", "MaAsLin2", "metagenomeSeq", "DACOMP", "ZicoSeq", "ANCOM",
             "ANCOMBC", "LinDA", "CAMP")
FDRs <- c(adapt_performance$FDR, aldex_performance$FDR, maaslin2_performance$FDR,
          metagenomeseq_performance$FDR, dacomp_performance$FDR, zicoseq_performance$FDR, 
          ancom_performance$FDR, ancombc_performance$FDR, linda_performance$FDR, camp_performance$FDR)
Powers <- c(adapt_performance$Power, aldex_performance$Power, maaslin2_performance$Power,
            metagenomeseq_performance$Power, dacomp_performance$Power, zicoseq_performance$Power, 
            ancom_performance$Power, ancombc_performance$Power, linda_performance$Power, camp_performance$Power)
Durations <- c(adapt_duration_noboot, aldex_duration, maaslin2_duration,
               metagenomeseq_duration, dacomp_duration, zicoseq_duration, ancom_duration, 
               ancombc_duration, linda_duration, camp_duration)
performance_summary <- data.frame(ID = myseed, Method=Methods,
                                  FDR = FDRs, Power=Powers, Duration=Durations)

output_filename <- sprintf("experiment_%d_%d.rds", choice, myseed)
saveRDS(performance_summary, file=file.path(folder, 'FoldChange', 'experiments', output_filename))



