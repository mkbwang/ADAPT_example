library(dada2); packageVersion("dada2")
library(ggplot2)
library(dplyr)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/16S/seqfetch"
# read metadata of the fastq files
metadata <- read.csv(file.path(folder, 'runinfo_Amplicon16S.csv'))
batchnumber <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

metadata_subset <- metadata %>% filter(batch == batchnumber)


# File parsing 
filtpathF <- file.path(folder, "filter_sequences", 'FWD') # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- file.path(folder, "filter_sequences", 'REV')
filtFs <- file.path(filtpathF, sprintf('%s_1.fastq.gz', metadata_subset$Run))
filtRs <- file.path(filtpathR, sprintf('%s_2.fastq.gz', metadata_subset$Run))
sample.names <- sapply(strsplit(basename(filtFs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names



# Learn forward error rates
set.seed(1)
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
#Plot & save errors 
Run_errF<-plotErrors(errF, nominalQ=TRUE)
Run_errR<-plotErrors(errR, nominalQ=TRUE)
ggsave(filename=paste(filtpathF, sprintf("run_errorFwd_%d.pdf", batchnumber), sep="/"), plot=Run_errF)
ggsave(filename=paste(filtpathR, sprintf("run_errorRev_%d.pdf", batchnumber), sep="/"), plot=Run_errR)


# denoise sequences and merge FWD and REV strands
mergers <- vector("list", length(sample.names))
dadaFs<-vector("list", length(sample.names))
dadaRs<-vector("list", length(sample.names))
names(mergers) <- sample.names
names(dadaFs)<-sample.names
names(dadaRs)<-sample.names
for(sam in sample.names) {
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  dadaFs[[sam]]<-ddF
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  dadaRs[[sam]]<-ddR
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}


# Construct sequence table 
seqtab_run <- makeSequenceTable(mergers)


saveRDS(seqtab_run, file.path(folder, "filter_sequences", sprintf("seqtable_%d.rds", batchnumber)))

