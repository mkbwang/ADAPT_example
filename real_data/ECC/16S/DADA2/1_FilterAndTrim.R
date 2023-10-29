library(dada2); packageVersion("dada2")

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/16S/seqfetch"
# File parse the raw sequences
pathF <- file.path(folder, "raw_sequences", "FWD")
pathR <- file.path(folder, "raw_sequences", "REV")
filtpathF <- file.path(folder, "filter_sequences", "FWD")
filtpathR <- file.path(folder, "filter_sequences", "REV")
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
subset <- seq((ID-1)*10+1, min(ID*10, length(fastqFs)))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")


RunOut<-filterAndTrim(fwd=file.path(pathF, fastqFs[subset]), filt=file.path(filtpathF, fastqFs[subset]),
              rev=file.path(pathR, fastqRs[subset]), filt.rev=file.path(filtpathR, fastqRs[subset]),
              trimLeft = c(8, 9), truncLen=c(240,210), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

save(RunOut, file=file.path(folder, "filter_sequences", sprintf("filter_trim_dada2_%d.Rdata", ID)))

