################################################
################################################SETUP
library(phyloseq)
library(stringr)
library(plyr)
library(dplyr)
library(decontam)
library(stringr)
library(ggplot2)
library(cowplot)

rm(list=ls())
folder <- '/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC'

data_folder <- file.path(folder, '16S', 'phyloseq')
metadata_folder <- file.path(folder, 'metadata')

# load the COHRA sample ID (private data)
load(file.path(metadata_folder, "SeqToCOHRAKey.Rdata")) # load metadata

# load the SRR information available on NCBI
amp_SRR_df <- read.csv(file.path(folder, '16S', 'seqfetch', 'runinfo_Amplicon16S.csv'))

## Match the sample name with the SRR ID, limit to true samples (no controls)
metadata_truesample <- metadata %>% filter(SampleType == "True Sample") %>%
  select(COHRAID, FoxCavID, DNAQuant, DNAQuant1000, AmplificationStatus) %>%
  dplyr::rename(Library.Name=FoxCavID)

## add an indicator for control samples
iscontrol <- grepl("control", amp_SRR_df$Library.Name)
amp_SRR_df <- amp_SRR_df %>% left_join(metadata_truesample, by="Library.Name")
amp_SRR_df$iscontrol <- grepl("control", amp_SRR_df$Library.Name)
amp_SRR_df$DNAQuant1000[is.na(amp_SRR_df$DNAQuant1000)] <- 1
amp_SRR_df$DNAQuant[is.na(amp_SRR_df$DNAQuant)] <- 0
rownames(amp_SRR_df) <- amp_SRR_df$Run
amp_SRR_df$Run <- NULL

# create phyloseq object
## count table name and the taxonomy name
taxname<-"tax_toSp_homd.rds"
seqname<-"seqtab_homd.rds"

seqtable <- readRDS(file.path(data_folder, seqname))
taxtable <- readRDS(file.path(data_folder, taxname))
## generate phyloseq object
phy.asv<- phyloseq(otu_table(seqtable, taxa_are_rows=FALSE),
                   tax_table(taxtable), sample_data(amp_SRR_df))
dna <- Biostrings::DNAStringSet(taxa_names(phy.asv))
names(dna) <- taxa_names(phy.asv)
phy.asv <- merge_phyloseq(phy.asv, dna)
taxa_names(phy.asv) <- paste0("ASV", seq(ntaxa(phy.asv)))
phy.asv <- subset_samples(phy.asv, sample_sums(phy.asv)>0)# filter out samples with zero reads
saveRDS(phy.asv, file=file.path(data_folder, "phyasv.rds"))

phy.asv <- readRDS(file.path(data_folder, "phyasv.rds"))

# decontamination
## set variables
decontam.method<-'either' #see ?isContaminant for method options
decontam.thresh<-c(0.1, 0.3) #see ?isContaminant for threshold

## decide if ASV is contaminant
contamdf<-isContaminant(phy.asv, batch="batch", method=decontam.method, neg="iscontrol", conc="DNAQuant1000", threshold=decontam.thresh)
phy.asv.nocontam <- subset_taxa(phy.asv, contamdf$contaminant==FALSE)


# filter samples and ASVs
readlimit<-1000

## a sample should be a true sample and have at least 1000 reads
phy.asv.truesample <- subset_samples(phy.asv.nocontam, !iscontrol & sample_sums(phy.asv.nocontam) > readlimit)
saveRDS(phy.asv.nocontam, file.path(data_folder, "phyasv_decontam.rds"))

## filter ASV
abund.filter<-0.05 # what abundance filter? [0 = 0% abundance to 1 = 100% abundance]
prev<-0.05 # what prevalence filter? [0=ASV present in 0% of samples to 1=ASV present in 100% of samples]
count_table <- otu_table(phy.asv.truesample)
count_mat <- count_table@.Data
relative_abundance <- count_mat / rowSums(count_mat)
maxRA <- apply(relative_abundance, 2, max)
prevalences <- colMeans(count_mat > 0)
taxa_filter <- maxRA > abund.filter | prevalences > prev
phy.asv.filter <- subset_taxa(phy.asv.truesample, taxa_filter)
metadata <- sample_data(phy.asv.filter)
sample_names(phy.asv.filter) <- metadata$COHRAID
saveRDS(phy.asv.filter, file.path(data_folder, "phyasv_decontam_filter.rds"))


