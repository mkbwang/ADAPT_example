
rm(list=ls())

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation"
suppressMessages(source(file.path(folder, "SparseDOSSA_setting.R")))

baseline_params$depth_fit[1] <- log(2e4)
baseline_params$depth_fit[2] <- 1


settings_df <- expand.grid(nSample=c(50, 100),
                           nTaxa=500,
                           balance_depth = c(FALSE, TRUE),
                           adj_depth = c(FALSE, TRUE))


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

sample_size <- settings_df$nSample[choice]
balance_depth <- settings_df$balance_depth[choice]
adj_depth <- settings_df$adj_depth[choice]

# choice equals 1 means similar sequencing depths at 2e4
# choice equals 2 means sequencing depths differ (1e4 and 2e4), and not adjusting for log sequencing depths
# choice equals 3 means sequencing depths differ (1e4 and 2e4), and adjusting for log sequencing depths


set.seed(myseed)
if (balance_depth){
  simulated_data <- SparseDOSSA2(template=baseline_params, n_sample=sample_size, n_feature=500,
                                 spike_metadata="none", median_read_depth=2e4)
  count_table <- simulated_data$simulated_data
} else{
  
  simulated_data <- SparseDOSSA2(template=baseline_params, n_sample=sample_size, n_feature=500,
                                    spike_metadata="none", median_read_depth=2e4)
  compositions <- simulated_data$simulated_matrices$rel
  count_table <- matrix(0, nrow=500, ncol=sample_size)
  seqdepths <- exp(c(rnorm(n=sample_size/2, mean=log(1e4), sd=1), rnorm(n=sample_size/2, mean=log(2e4), sd=1))) |> 
    round()
  for (j in 1:sample_size){
    count_table[,j] <- rmultinom(n=1, size=seqdepths[j], prob=compositions[, j])
  }
}

metadata_df <- data.frame(main=rep(c(0, 1), each=sample_size/2), confounder=rnorm(n=sample_size, mean=0, sd=1))
colnames(count_table) <- sprintf("Sample%d", seq(1, sample_size))
rownames(count_table) <- sprintf("Taxon%d", seq(1, nrow(count_table)))
rownames(metadata_df) <- sprintf("Sample%d", seq(1, sample_size))


prevalence_cutoff <- 0.05
seqdepth_cutoff <- 1e3
prevalences <- rowMeans(count_table > 0)
taxa_filter <- prevalences > 0.05
count_table <- count_table[taxa_filter, ]


seqdepths <- colSums(count_table)


seqdepth_filter <- seqdepths > seqdepth_cutoff
count_table <- count_table[, seqdepth_filter, drop=F]
metadata_df <- metadata_df[seqdepth_filter, ]
metadata_df$logdepth <- log(colSums(count_table))



# ADAPT
suppressMessages(library(ADAPT))
begin <- proc.time()
if (adj_depth){
  adapt_output <- adapt(otu_table=count_table, metadata=metadata_df,
                        covar="main", adjust=c("logdepth"), prevalence_cutoff=0.05, depth_cutoff=0,
                        taxa_are_rows = T, boot=F)
} else{
  adapt_output <- adapt(otu_table=count_table, metadata=metadata_df,
                        covar="main", adjust=NULL, prevalence_cutoff=0.05, depth_cutoff=0,
                        taxa_are_rows = T, boot=F)
}

adapt_time <- proc.time() - begin
adapt_duration_noboot <- adapt_time[3]
source(file.path(folder, "methods", "adapt_utils.R"))
adapt_performance_noboot <- suppressMessages(evaluation_adapt(adapt_result=adapt_output,
                                                       nullcase=T))



# Aldex2
suppressMessages(library(ALDEx2))
begin <- proc.time()
if (adj_depth){
  mm <- model.matrix(~main+logdepth, metadata_df)
  aldex_perm <- suppressMessages(aldex.clr(count_table, mm, mc.samples=4, denom="all"))
  aldex_result <- aldex.glm(aldex_perm, mm)
} else{
  aldex_result <- suppressMessages(aldex(reads=count_table, conditions=metadata_df$main, mc.samples=128,
                                         test="t"))
}
aldex_time <- proc.time() - begin
aldex_duration <- aldex_time[3]
source(file.path(folder, "methods", "aldex2_utils.R"))
if (adj_depth){
  aldex_performance <- suppressMessages(evaluation_aldex2(aldex_result=aldex_result, test="glm",
                                                          nullcase = T))
} else{
  aldex_performance <- suppressMessages(evaluation_aldex2(aldex_result=aldex_result, test="ttest",
                                                          nullcase = T))
}


# Maaslin2
suppressMessages(library(Maaslin2))
count_df <- data.frame(t(count_table))
begin <- proc.time()
if (adj_depth){
  maaslin2_output <-invisible(Maaslin2(input_data=count_df, input_metadata=metadata_df,
                                       output = file.path(folder, 'Null_Case', 'maaslin2_output', sprintf("seed_%d", myseed)),
                                       min_prevalence=0.05,
                                       max_significance = 0.05,
                                       fixed_effects=c("main", "logdepth"),
                                       cores=4, plot_heatmap = F,plot_scatter = F))
} else{
  maaslin2_output <-invisible(Maaslin2(input_data=count_df, input_metadata=metadata_df,
                                       output = file.path(folder, 'Null_Case', 'maaslin2_output', sprintf("seed_%d", myseed)),
                                       min_prevalence=0.05,
                                       max_significance = 0.05,
                                       fixed_effects=c("main"),
                                       cores=4, plot_heatmap = F,plot_scatter = F))
}
maaslin2_time <- proc.time() - begin
maaslin2_duration <- maaslin2_time[3]
maaslin2_result <- maaslin2_output$results
source(file.path(folder, "methods", "maaslin2_utils.R"))
maaslin2_performance <- suppressMessages(evaluation_maaslin2(maaslin2_result=maaslin2_result,
                                                             nullcase=T))




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
source(file.path(folder, "methods", "metagenomeseq_utils.R"))
metagenomeseq_performance <- suppressMessages(evaluation_metagenomeseq(MR_result=MRfit, nullcase=T))



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
source(file.path(folder, "methods", "DACOMP_utils.R"))
dacomp_performance <- evaluation_dacomp(dacomp_result=dacomp_output,
                                        nullcase=T)



# ZicoSeq
library(GUniFrac)
begin <- proc.time()
if (adj_depth){
  zicoseq_output <- ZicoSeq(meta.dat=metadata_df,
                            feature.dat=count_table,
                            grp.name="main", adj.name="logdepth",
                            is.fwer=T,
                            perm.no=999)
} else{
  zicoseq_output <- ZicoSeq(meta.dat=metadata_df,
                            feature.dat=count_table,
                            grp.name="main",
                            is.fwer=T,
                            perm.no=999)
}
zicoseq_result <- cbind.data.frame(p.raw = zicoseq_output$p.raw,
                                   p.adj.fdr = zicoseq_output$p.adj.fdr,
                                   p.adj.fwer = zicoseq_output$p.adj.fwer)
zicoseq_time <- proc.time()-begin
zicoseq_duration <- zicoseq_time[3]
source(file.path(folder, "methods", "zicoseq_utils.R"))
zicoseq_performance <- evaluation_zicoseq(zicoseq_result=zicoseq_result,
                                          nullcase=T)



# RDB
# library(RDB)
# 
# begin <- proc.time()
# depths <- rowSums(count_table)
# compositions <- count_table / depths
# rdb_output <- rdb(P=compositions, Z=metadata$main,
#                   alpha=0.05, fdr=F)
# rdb_time <- proc.time() - begin
# rdb_duration <- rdb_time[3]
# source(file.path(folder, "methods", "RDB_utils.R"))
# rdb_performance <- evaluation_rdb(taxa_truth=taxa_info,
#                                   rdb_result=rdb_output,
#                                   nullcase = T)



# LOCOM
# library(LOCOM)
# ptm <- proc.time()
# locom_output <- locom(otu.table=count_table,
#                       Y=metadata$main, C=metadata$confounder, fdr.nominal=0.05, prev.cut=0.05,
#                       n.perm.max=1e4, seed=1, Firth.thresh = 1,
#                       n.cores=4)
# locom_time <- proc.time() - ptm
# locom_duration <- locom_time[3]
# source(file.path(folder, "methods", "locom_utils.R"))
# locom_performance <- evaluation_locom(taxa_info, locom_result=locom_output,
#                                       nullcase=T)


# ANCOMBC
suppressMessages(library(ANCOMBC))
suppressMessages(library(phyloseq))
phyobj <- phyloseq(otu_table(count_table, taxa_are_rows = T),
                   sample_data(metadata_df))
ptm <- proc.time()
if (adj_depth){
  ancombc_output <- ancombc2(phyobj, fix_formula = 'main+logdepth', p_adj_method='BH', global=T,
                             group='main', struc_zero=FALSE, prv_cut=0.05, n_cl=4)  
} else{
  ancombc_output <- ancombc2(phyobj, fix_formula = 'main', p_adj_method='BH', global=T,
                             group='main', struc_zero=FALSE, prv_cut=0.05, n_cl=4)  
}
ancombc_time <- proc.time() - ptm
ancombc_duration <- ancombc_time[3]
source(file.path(folder, "methods", "ancombc_utils.R"))
ancombc_performance <- evaluation_ancombc(ancombc_result=ancombc_output,
                                          nullcase=T)


# LinDA
suppressMessages(library(MicrobiomeStat))
begin <- proc.time()
if (adj_depth){
  linda_output <- linda(feature.dat=count_table,
                        meta.dat=metadata_df,
                        formula="~main+logdepth",
                        n.cores=4)
} else{
  linda_output <- linda(feature.dat=count_table,
                        meta.dat=metadata_df,
                        formula="~main",
                        n.cores=4)
}
linda_time <- proc.time() - begin
linda_duration <- linda_time[3]
source(file.path(folder, "methods", "linda_utils.R"))
linda_performance <- evaluation_linda(linda_result=linda_output,
                                      nullcase=T)



# RioNorm2
# library(RioNorm2)
# 
# begin <- proc.time()
# rionorm2_output <-RioNorm2_test(OTU_table=t(count_table),
#                       class=metadata$main,FDR=1) |> as.data.frame()
# rionorm2_output$padj <- as.numeric(rionorm2_output$padj)
# rionorm2_output$pvalue <- as.numeric(rionorm2_output$pvalue)
# rionorm2_time <- proc.time() - begin
# rionorm2_duration <- rionorm2_time[3]
# source(file.path(folder, "methods", "RioNorm2_utils.R"))
# rionorm2_performance <- evaluation_rionorm2(taxa_truth=taxa_info,
#                                       rionorm2_result=rionorm2_output,
#                                       nullcase=F)



# aggregate all the performances
Methods <- c("ADAPT", "ALDEx2", "Maaslin2", "metagenomeSeq", "DACOMP", "ZicoSeq",
             "ANCOMBC", "LinDA")
FPR_1s <- c(adapt_performance_noboot$FPR_1, aldex_performance$FPR_1, maaslin2_performance$FPR_1,
          metagenomeseq_performance$FPR_1, dacomp_performance$FPR_1, zicoseq_performance$FPR_1, 
          ancombc_performance$FPR_1, linda_performance$FPR_1)
FPR_5s <- c(adapt_performance_noboot$FPR_5, aldex_performance$FPR_5, maaslin2_performance$FPR_5,
            metagenomeseq_performance$FPR_5, dacomp_performance$FPR_5, zicoseq_performance$FPR_5, 
            ancombc_performance$FPR_5, linda_performance$FPR_5)
FPR_10s <- c(adapt_performance_noboot$FPR_10, aldex_performance$FPR_10, maaslin2_performance$FPR_10,
            metagenomeseq_performance$FPR_10, dacomp_performance$FPR_10, zicoseq_performance$FPR_10, 
            ancombc_performance$FPR_10, linda_performance$FPR_10)
FWE_1s <- c(adapt_performance_noboot$FWE_1, aldex_performance$FWE_1, maaslin2_performance$FWE_1,
            metagenomeseq_performance$FWE_1, dacomp_performance$FWE_1, zicoseq_performance$FWE_1, 
            ancombc_performance$FWE_1, linda_performance$FWE_1)
FWE_5s <- c(adapt_performance_noboot$FWE_5, aldex_performance$FWE_5, maaslin2_performance$FWE_5,
            metagenomeseq_performance$FWE_5, dacomp_performance$FWE_5, zicoseq_performance$FWE_5, 
            ancombc_performance$FWE_5, linda_performance$FWE_5)
FWE_10s <- c(adapt_performance_noboot$FWE_10, aldex_performance$FWE_10, maaslin2_performance$FWE_10,
             metagenomeseq_performance$FWE_10, dacomp_performance$FWE_10, zicoseq_performance$FWE_10, 
             ancombc_performance$FWE_10, linda_performance$FWE_10)
Durations <- c(adapt_duration_noboot, aldex_duration, maaslin2_duration,
               metagenomeseq_duration, dacomp_duration, zicoseq_duration, ancombc_duration, linda_duration)
performance_summary <- data.frame(ID = myseed, Method=Methods,
                                  FPR_1=FPR_1s, FPR_5=FPR_5s, FPR_10=FPR_10s,
                                  FWE_1=FWE_1s, FWE_5=FWE_5s, FWE_10=FWE_10s,
                                  Duration=Durations)

output_filename <- sprintf("experiment_%d_%d.rds", choice, myseed)
saveRDS(performance_summary, file=file.path(folder, 'Null_Case', 'experiments', output_filename))

