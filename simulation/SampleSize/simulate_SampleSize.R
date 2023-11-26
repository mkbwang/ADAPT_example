
rm(list=ls())

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/simulation"
source(file.path(folder, "simulation_utils.R"))
dirparams <- readRDS(file.path(folder, "DM_param.rds"))

settings_df <- expand.grid(nSample=c(50, 100, 200),
                           seqdepth_mean=2e4,
                           nTaxa=1000,
                           depth_confound=FALSE,
                           propDA=0.1,
                           main_mean_logfold=log(4),
                           direction=c("balanced", "unbalanced"),
                           propcf=0,
                           cf_sd=1,
                           cf_main_corr=0,
                           cf_mean_logfold=1,
                           zero_inflation=0,
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


simulated_data <- SimulateCount(dirparams, 
                                nSample = settings_df$nSample[choice],
                                nTaxa = settings_df$nTaxa[choice],
                                seqdepth_mean=settings_df$seqdepth_mean[choice], 
                                seqdepth_theta=10, depth_confound=settings_df$depth_confound[choice],
                                propDA=settings_df$propDA[choice], 
                                main_mean_logfold = settings_df$main_mean_logfold[choice], 
                                main_sd_logfold=0.1, 
                                direction = settings_df$direction[choice],
                                propcf = settings_df$propcf[choice], 
                                cf_sd = settings_df$cf_sd[choice], 
                                cf_main_corr = settings_df$cf_main_corr[choice], 
                                cf_mean_logfold = settings_df$cf_mean_logfold[choice], 
                                cf_sd_logfold=0.1,
                                zero_inflation=settings_df$zero_inflation[choice], 
                                seed=myseed)

count_table <- simulated_data$count_mat
metadata <- simulated_data$sample_metadata
taxa_info <- simulated_data$taxa_info

prevalence_cutoff <- 0.05
prevalences <- colMeans(count_table > 0)
taxa_filter <- prevalences > 0.05
count_table <- count_table[, taxa_filter]
taxa_info <- taxa_info[taxa_filter, , drop=F]




# ADAPT
suppressMessages(library(ADAPT))
begin <- proc.time()
adapt_output <- adapt(otu_table=count_table, metadata=metadata,
                      covar="main",  prevalence_cutoff=0, depth_cutoff=0,
                      taxa_are_rows = F, boot=F)
adapt_time <- proc.time() - begin
adapt_duration_noboot <- adapt_time[3]
source(file.path(folder, "methods", "adapt_utils.R"))
adapt_performance_noboot <- suppressMessages(evaluation_adapt(taxa_truth=taxa_info,
                                                       adapt_result=adapt_output,
                                                       nullcase=F))

begin <- proc.time()
adapt_output <- adapt(otu_table=count_table, metadata=metadata,
                      covar="main", adjust="confounder", prevalence_cutoff=0.05, depth_cutoff=0,
                      taxa_are_rows = F, boot=T)
adapt_time <- proc.time() - begin
adapt_duration_boot <- adapt_time[3]
adapt_performance_boot <- suppressMessages(evaluation_adapt(taxa_truth=taxa_info,
                                                       adapt_result=adapt_output,
                                                       nullcase=F))



# Aldex2
suppressMessages(library(ALDEx2))
begin <- proc.time()
mm <- model.matrix(~main+confounder, metadata)
aldex_perm <- suppressMessages(aldex.clr(t(count_table), mm, mc.samples=4, denom="all"))
aldex_result <- aldex.glm(aldex_perm, mm, fdr.method="BH")
aldex_time <- proc.time() - begin
aldex_duration <- aldex_time[3]
source(file.path(folder, "methods", "aldex2_utils.R"))
aldex_performance <- suppressMessages(evaluation_aldex2(taxa_truth=taxa_info, 
                                       aldex_result=aldex_result,
                                       nullcase = F))


# MetagenomeSeq
suppressMessages(library(metagenomeSeq))
MRmetadata <- AnnotatedDataFrame(metadata)
MRobj <- newMRexperiment(counts=t(count_table),
                         phenoData=MRmetadata)
begin <- proc.time()
MRobj <- wrenchNorm(MRobj, condition=metadata$main)
mod <- model.matrix(~main, data=pData(MRobj))
MRfit <- fitFeatureModel(MRobj, mod)
metagenomeseq_time <- proc.time() - begin
metagenomeseq_duration <- metagenomeseq_time[3]
source(file.path(folder, "methods", "metagenomeseq_utils.R"))
metagenomeseq_performance <- suppressMessages(evaluation_metagenomeseq(taxa_truth=taxa_info,
                                           MR_result=MRfit, nullcase=F))



# DACOMP
library(dacomp)
begin <- proc.time()
selected_references <- dacomp.select_references(
  X = count_table, verbose = F, run_in_parallel = T, Nr.Cores=4,
  minimal_TA=10)
updated_selected_references <- dacomp.validate_references(X = count_table, #Counts
                                                          Y = metadata$main, #Traits
                                                          ref_obj = selected_references, #reference checked, must be object from dacomp.select_references(...)
                                                          test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                          Q_validation = 0.1, #FDR level for checking
                                                          Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                          Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                          NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                          Verbose = F)
dacomp_output <- dacomp.test(X = count_table, #counts data
                          y = metadata$main, #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = updated_selected_references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          disable_DSFDR = T,
                          q=0.05,
                          verbose = F)

dacomp_time <- proc.time() - begin
dacomp_duration <- dacomp_time[3]
source(file.path(folder, "methods", "DACOMP_utils.R"))
dacomp_performance <- evaluation_dacomp(taxa_truth=taxa_info,
                                        dacomp_result=dacomp_output,
                                        nullcase=F)



# ZicoSeq
library(GUniFrac)
begin <- proc.time()
zicoseq_output <- ZicoSeq(meta.dat=metadata,
                          feature.dat=t(count_table),
                          grp.name="main",
                          adj.name="confounder",is.fwer=T,
                          perm.no=999)
zicoseq_result <- cbind.data.frame(p.raw = zicoseq_output$p.raw,
                                p.adj.fdr = zicoseq_output$p.adj.fdr,
                                p.adj.fwer = zicoseq_output$p.adj.fwer)
zicoseq_time <- proc.time()-begin
zicoseq_duration <- zicoseq_time[3]
source(file.path(folder, "methods", "zicoseq_utils.R"))
zicoseq_performance <- evaluation_zicoseq(taxa_truth=taxa_info,
                                          zicoseq_result=zicoseq_result,
                                          nullcase=F)



# RDB
library(RDB)

begin <- proc.time()
depths <- rowSums(count_table)
compositions <- count_table / depths
rdb_output <- rdb(P=compositions, Z=metadata$main,
                  alpha=0.05, fdr=T)
rdb_time <- proc.time() - begin
rdb_duration <- rdb_time[3]
source(file.path(folder, "methods", "RDB_utils.R"))
rdb_performance <- evaluation_rdb(taxa_truth=taxa_info,
                                  rdb_result=rdb_output,
                                  nullcase = F)



# LOCOM
library(LOCOM)
ptm <- proc.time()
locom_output <- locom(otu.table=count_table,
                      Y=metadata$main, C=metadata$confounder, fdr.nominal=0.05, prev.cut=0.05,
                      n.perm.max=1e4, seed=1, Firth.thresh = 1,
                      n.cores=4)
locom_time <- proc.time() - ptm
locom_duration <- locom_time[3]
source(file.path(folder, "methods", "locom_utils.R"))
locom_performance <- evaluation_locom(taxa_info, locom_result=locom_output,
                                      nullcase=F)


# ANCOMBC
suppressMessages(library(ANCOMBC))
suppressMessages(library(phyloseq))
phyobj <- phyloseq(otu_table(count_table, taxa_are_rows = F),
                   sample_data(metadata))
ptm <- proc.time()
ancombc_output <- ancombc2(phyobj, fix_formula = 'main + confounder', p_adj_method='BH', global=T,
                          group='main', struc_zero=FALSE, prv_cut=0.05, n_cl=4)
ancombc_time <- proc.time() - ptm
ancombc_duration <- ancombc_time[3]
source(file.path(folder, "methods", "ancombc_utils.R"))
ancombc_performance <- evaluation_ancombc(taxa_truth=taxa_info,
                                          ancombc_result=ancombc_output,
                                          nullcase=F)


# LinDA
suppressMessages(library(LinDA))

begin <- proc.time()
linda_output <- linda(otu.tab=t(count_table),
                      meta=metadata,
                      formula="~main+confounder",
                      n.cores=4)
linda_time <- proc.time() - begin
linda_duration <- linda_time[3]
source(file.path(folder, "methods", "linda_utils.R"))
linda_performance <- evaluation_linda(taxa_truth=taxa_info,
                                      linda_result=linda_output,
                                      nullcase=F)


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
Methods <- c("ADAPT_noboot", "ADAPT_boot", "ALDEx2", "MetagenomeSeq", "DACOMP", "ZicoSeq", "RDB", "LOCOM", "ANCOMBC", "LinDA")
FDRs <- c(adapt_performance_noboot$FDR, adapt_performance_boot$FDR, aldex_performance$FDR, metagenomeseq_performance$FDR,
         dacomp_performance$FDR, zicoseq_performance$FDR, rdb_performance$FDR, locom_performance$FDR,
         ancombc_performance$FDR, linda_performance$FDR)
FPRs <- c(adapt_performance_noboot$FPR, adapt_performance_boot$FPR, aldex_performance$FPR, metagenomeseq_performance$FPR,
          dacomp_performance$FPR, zicoseq_performance$FPR, rdb_performance$FPR, locom_performance$FPR,
          ancombc_performance$FPR, linda_performance$FPR)
Powers <- c(adapt_performance_noboot$Power, adapt_performance_boot$Power, aldex_performance$Power, metagenomeseq_performance$Power,
            dacomp_performance$Power, zicoseq_performance$Power, rdb_performance$Power, locom_performance$Power,
            ancombc_performance$Power, linda_performance$Power)
Durations <- c(adapt_duration_noboot, adapt_duration_boot, aldex_duration, metagenomeseq_duration, dacomp_duration, zicoseq_duration,
              rdb_duration, locom_duration, ancombc_duration, linda_duration)
performance_summary <- data.frame(ID = myseed, Method=Methods,
                                  FDR = FDRs,FPR=FPRs, Power=Powers, Duration=Durations)

output_filename <- sprintf("experiment_%d_%d.rds", choice, myseed)
saveRDS(performance_summary, file=file.path(folder, 'SampleSize', 'experiments', output_filename))

