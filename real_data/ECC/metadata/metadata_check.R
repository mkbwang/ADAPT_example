rm(list=ls())
library(dplyr)

ECC_folder <- "real_data/ECC"

# https://www.ncbi.nlm.nih.gov/Traces/study/ PRJNA752888
SRA_metadata <- read.table(file.path(ECC_folder, "metadata","SraRunTable.txt"), header=TRUE,
                           sep=",")
metadata_subset <- SRA_metadata %>% select(Run, Assay.Type, Library.Name, dbGaP_ID,
                                           batch, Host_age, host_sex, geo_loc_name, env_medium,
                                           Collection_Date)

WGS_metadata <- metadata_subset %>% filter(Assay.Type == "WGS")
Amplicon16S_metadata <- metadata_subset %>% filter(Assay.Type == "AMPLICON")

write.csv(WGS_metadata, file.path(ECC_folder, 'WGS', 'seqfetch', 'runinfo_WGS.csv'),
          row.names=FALSE, quote=FALSE)

write.csv(Amplicon16S_metadata, file.path(ECC_folder, '16S', 'seqfetch', 'runinfo_Amplicon16S.csv'),
          row.names=FALSE, quote=FALSE)

