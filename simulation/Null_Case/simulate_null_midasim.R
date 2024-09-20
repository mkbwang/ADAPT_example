
rm(list=ls())

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation"
suppressMessages(source(file.path(folder, "SparseDOSSA_setting.R")))
suppressMessages(source(file.path(folder, "MIDASim_setting.R")))
suppressMessages(source(file.path(folder, "camp.R")))

baseline_params <- readRDS(file.path(folder, 'midasim_baseline.rds'))

baseline_params$depth_fit[1] <- log(1e4)
baseline_params$depth_fit[2] <- 1


settings_df <- expand.grid(nSample=c(50, 100),
                           nTaxa=500,
                           depth_fold=c(1,4,10))


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
taxa_num <- 500
depth_fold <- settings_df$depth_fold[choice]



set.seed(myseed)
counts_group1 <- midasim_null_simulate(baseline=baseline_params,
                                       num_sample=sample_size/2,
                                       num_taxa=taxa_num)
baseline_params$depth_fit[1] <- baseline_params$depth_fit[1] + log(depth_fold)
counts_group2 <- midasim_null_simulate(baseline=baseline_params,
                                       num_sample=sample_size/2,
                                       num_taxa=taxa_num)
count_table <- rbind(counts_group1, counts_group2)


metadata_df <- data.frame(main=rep(c(0, 1), each=sample_size/2), confounder=rnorm(n=sample_size, mean=0, sd=1))
rownames(count_table) <- sprintf("Sample%d", seq(1, sample_size))
colnames(count_table) <- sprintf("Taxon%d", seq(1, taxa_num))
rownames(metadata_df) <- sprintf("Sample%d", seq(1, sample_size))


prevalence_cutoff <- 0.05

seqdepth_cutoff <- 1e3
prevalences <- colMeans(count_table > 0)
taxa_filter <- prevalences > 0.05
count_table <- count_table[, taxa_filter]


seqdepths <- rowSums(count_table)


seqdepth_filter <- seqdepths > seqdepth_cutoff
count_table <- count_table[seqdepth_filter, ,drop=F]
metadata_df <- metadata_df[seqdepth_filter, ]


suppressMessages(library(phyloseq))
phyobj <- phyloseq(otu_table(count_table, taxa_are_rows = F),
                   sample_data(metadata_df))


# ADAPT
suppressMessages(library(ADAPT))
begin <- proc.time()
adapt_output <- adapt(input_data=phyobj, cond.var="main")
adapt_time <- proc.time() - begin
adapt_duration_noboot <- adapt_time[3]
source(file.path(folder, "methods_evaluation", "adapt_utils.R"))
adapt_performance <- suppressMessages(evaluation_adapt(adapt_result=adapt_output,
                                                       nullcase=T))



# Aldex2
suppressMessages(library(ALDEx2))
begin <- proc.time()
aldex_result <- suppressMessages(aldex(reads=t(count_table), conditions=metadata_df$main, mc.samples=128,
                                         test="t"))

aldex_time <- proc.time() - begin
aldex_duration <- aldex_time[3]
source(file.path(folder, "methods_evaluation", "aldex2_utils.R"))
aldex_performance <- suppressMessages(evaluation_aldex2(aldex_result=aldex_result, test="ttest",
                                                        nullcase = T))


# Maaslin2
suppressMessages(library(Maaslin2))
count_df <- data.frame(count_table)
begin <- proc.time()
maaslin2_output <-invisible(Maaslin2(input_data=count_df, input_metadata=metadata_df,
                                     output = file.path(folder, 'Null_Case', 'maaslin2_output', sprintf("seed_%d", myseed)),
                                     min_prevalence=0.05,
                                     max_significance = 0.05,
                                     fixed_effects=c("main"),
                                     cores=4, plot_heatmap = F,plot_scatter = F))
maaslin2_time <- proc.time() - begin
maaslin2_duration <- maaslin2_time[3]
maaslin2_result <- maaslin2_output$results
source(file.path(folder, "methods_evaluation", "maaslin2_utils.R"))
maaslin2_performance <- suppressMessages(evaluation_maaslin2(maaslin2_result=maaslin2_result,
                                                             nullcase=T))




# MetagenomeSeq
suppressMessages(library(metagenomeSeq))
MRmetadata <- AnnotatedDataFrame(metadata_df)
MRobj <- newMRexperiment(counts=t(count_table),
                         phenoData=MRmetadata)
begin <- proc.time()
MRobj <- wrenchNorm(MRobj, condition=metadata_df$main)
mod <- model.matrix(~main, data=pData(MRobj))
MRfit <- fitFeatureModel(MRobj, mod)
metagenomeseq_time <- proc.time() - begin
metagenomeseq_duration <- metagenomeseq_time[3]
source(file.path(folder, "methods_evaluation", "metagenomeseq_utils.R"))
metagenomeseq_performance <- suppressMessages(evaluation_metagenomeseq(MR_result=MRfit, nullcase=T))



# DACOMP
library(dacomp)
begin <- proc.time()
selected_references <- dacomp.select_references(
  X = count_table, verbose = F, run_in_parallel = T, Nr.Cores=4,
  minimal_TA=1)
updated_selected_references <- dacomp.validate_references(X = count_table, #Counts
                                                          Y = metadata_df$main, #Traits
                                                          ref_obj = selected_references, #reference checked, must be object from dacomp.select_references(...)
                                                          test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                          Q_validation = 0.1, #FDR level for checking
                                                          Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                          Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                          NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                          Verbose = F)
dacomp_output <- dacomp.test(X = count_table, #counts data
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
dacomp_performance <- evaluation_dacomp(dacomp_result=dacomp_output,
                                        nullcase=T)



# ZicoSeq
library(GUniFrac)
begin <- proc.time()
zicoseq_output <- ZicoSeq(meta.dat=metadata_df,
                          feature.dat=t(count_table),
                          grp.name="main",
                          is.fwer=T,
                          perm.no=999)
zicoseq_result <- cbind.data.frame(p.raw = zicoseq_output$p.raw,
                                   p.adj.fdr = zicoseq_output$p.adj.fdr,
                                   p.adj.fwer = zicoseq_output$p.adj.fwer)
zicoseq_time <- proc.time()-begin
zicoseq_duration <- zicoseq_time[3]
source(file.path(folder, "methods_evaluation", "zicoseq_utils.R"))
zicoseq_performance <- evaluation_zicoseq(zicoseq_result=zicoseq_result,
                                          nullcase=T)



# ANCOMBC
suppressMessages(library(ANCOMBC))
ptm <- proc.time()
ancombc_output <- ancombc2(phyobj, fix_formula = 'main', p_adj_method='BH', global=T,
                           group='main', struc_zero=FALSE, prv_cut=0.05, n_cl=4)  
ancombc_time <- proc.time() - ptm
ancombc_duration <- ancombc_time[3]
source(file.path(folder, "methods_evaluation", "ancombc_utils.R"))
ancombc_performance <- evaluation_ancombc(ancombc_result=ancombc_output,
                                          nullcase=T)


# LinDA
suppressMessages(library(MicrobiomeStat))
begin <- proc.time()
linda_output <- linda(feature.dat=t(count_table),
                      meta.dat=metadata_df,
                      formula="~main",
                      n.cores=4)
linda_time <- proc.time() - begin
linda_duration <- linda_time[3]
source(file.path(folder, "methods_evaluation", "linda_utils.R"))
linda_performance <- evaluation_linda(linda_result=linda_output,
                                      nullcase=T)

# CAMP
begin <- proc.time()
camp_output <- camp(X=count_table, cov=metadata_df)
camp_time <- proc.time() - begin
camp_duration <- camp_time[3]
source(file.path(folder, "methods_evaluation", "camp_utils.R"))
camp_performance <- evaluation_camp(camp_result = camp_output,
                                    nullcase=TRUE)

# aggregate all the performances
Methods <- c("ADAPT", "ALDEx2", "Maaslin2", "metagenomeSeq", "DACOMP", "ZicoSeq",
             "ANCOMBC", "LinDA", "CAMP")
FPR_1s <- c(adapt_performance$FPR_1, aldex_performance$FPR_1, maaslin2_performance$FPR_1,
          metagenomeseq_performance$FPR_1, dacomp_performance$FPR_1, zicoseq_performance$FPR_1, 
          ancombc_performance$FPR_1, linda_performance$FPR_1, camp_performance$FPR_1)
FPR_5s <- c(adapt_performance$FPR_5, aldex_performance$FPR_5, maaslin2_performance$FPR_5,
            metagenomeseq_performance$FPR_5, dacomp_performance$FPR_5, zicoseq_performance$FPR_5, 
            ancombc_performance$FPR_5, linda_performance$FPR_5, camp_performance$FPR_5)
FPR_10s <- c(adapt_performance$FPR_10, aldex_performance$FPR_10, maaslin2_performance$FPR_10,
            metagenomeseq_performance$FPR_10, dacomp_performance$FPR_10, zicoseq_performance$FPR_10, 
            ancombc_performance$FPR_10, linda_performance$FPR_10, camp_performance$FPR_10)
num_FDs <- c(adapt_performance$num_FD, aldex_performance$num_FD, maaslin2_performance$num_FD,
            metagenomeseq_performance$num_FD, dacomp_performance$num_FD, zicoseq_performance$num_FD, 
            ancombc_performance$num_FD, linda_performance$num_FD, camp_performance$num_FD)

Durations <- c(adapt_duration_noboot, aldex_duration, maaslin2_duration,
               metagenomeseq_duration, dacomp_duration, zicoseq_duration, ancombc_duration, linda_duration, camp_duration)
performance_summary <- data.frame(ID = myseed, Method=Methods,
                                  FPR_1=FPR_1s, FPR_5=FPR_5s, FPR_10=FPR_10s,
                                  num_FD=num_FDs, Duration=Durations)

output_filename <- sprintf("experiment_%d_%d.rds", choice, myseed)
saveRDS(performance_summary, file=file.path(folder, 'Null_Case', 'experiments_midasim', output_filename))

