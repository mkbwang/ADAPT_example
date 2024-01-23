rm(list=ls())

# Differential Abundance Analysis
library(dplyr)
folder <- 'real_data/ECC/'

# load the processed phyloseq data
ECC_16S <- readRDS(file.path(folder, "16S/phyloseq/phyasv_decontam_filter.rds"))
count_matrix <- otu_table(ECC_16S)@.Data
taxonomies <- tax_table(ECC_16S)@.Data
sample_metadata <- data.frame(sample_data(ECC_16S)) 
sample_metadata <- sample_metadata %>% dplyr::select(COHRAID, Library.Name, dbGaP_ID,
                                              host_sex, geo_loc_name,
                                              AmplificationStatus)

# external metadata
load(file.path(folder, "metadata", "MetaVisit.Rdata"))
external_metadata <- MetaVisit
external_columns <- colnames(external_metadata)
external_metadata <- external_metadata %>% dplyr::select(COHRAID, Visit, IncidentVisit, BabySubjectID, MotherSubjectID,
                                                  Delivery, BabySex, AgeAtExamMonths, CaseStatus, CaseEver)


# all visits
total_sample_metadata <- sample_metadata %>% select(-host_sex) %>% inner_join(external_metadata, by="COHRAID")
## two individuals have duplicate/conflicting metadata
intermediate1 <- total_sample_metadata %>% filter(BabySubjectID == 11000244) %>% filter(CaseEver == "Case")
intermediate2 <- total_sample_metadata %>% filter(BabySubjectID == 22000022) %>% filter(IncidentVisit == 10)
intermediate3 <- total_sample_metadata %>% filter(BabySubjectID != 11000244 & BabySubjectID != 22000022) 
total_sample_metadata_deduplicate <- rbind(intermediate1, intermediate2, intermediate3)

rownames(total_sample_metadata_deduplicate) <- total_sample_metadata_deduplicate$COHRAID
sample_data(ECC_16S) <- total_sample_metadata_deduplicate


saveRDS(ECC_16S, file="real_data/ECC/16S/phyloseq/ECC_16S_cleaned.rds")

# filter out visits at age 12 months (visit 5)
external_visit_12 <- external_metadata %>% filter(Visit == 5)

## merge the sample information from SRA and extra metadata
samples_visit12 <- sample_metadata %>% inner_join(external_visit_12, by="COHRAID")

## some individuals have conflicting information about the latest visit, deduplicate
samples_visit12_nocase <- samples_visit12 %>% filter(CaseStatus != "Case (incident visit)") 
latest_visit <- samples_visit12_nocase %>% 
  group_by(COHRAID) %>% 
  dplyr::summarise(latest = max(IncidentVisit))
samples_visit12_nocase <- samples_visit12_nocase %>% inner_join(latest_visit, by="COHRAID")
samples_visit12_nocase <- samples_visit12_nocase %>% filter(IncidentVisit == latest)
samples_visit12_nocase$latest <- NULL
rownames(samples_visit12_nocase) <- samples_visit12_nocase$COHRAID

ECC16S_month12 <- subset_samples(ECC_16S, COHRAID %in% samples_visit12_nocase$COHRAID)
sample_data(ECC16S_month12) <- samples_visit12_nocase

saveRDS(ECC16S_month12, file=file.path(folder, "16S/phyloseq/phyasv_visit12.rds"))


# filter out visits at age 24 months (visit 7)
external_visit_24 <- external_metadata %>% filter(Visit == 7)


## merge the sample information from SRA and extra metadata
samples_visit24 <- sample_metadata %>% inner_join(external_visit_24, by="COHRAID")

## some individuals have conflicting information about the latest visit, deduplicate
samples_visit24_nocase <- samples_visit24 %>% filter(CaseStatus != "Case (incident visit)") 
latest_visit <- samples_visit24_nocase %>% 
  group_by(COHRAID) %>% 
  dplyr::summarise(latest = max(IncidentVisit))
samples_visit24_nocase <- samples_visit24_nocase %>% inner_join(latest_visit, by="COHRAID")
samples_visit24_nocase <- samples_visit24_nocase %>% filter(IncidentVisit == latest)
samples_visit24_nocase$latest <- NULL
rownames(samples_visit24_nocase) <- samples_visit24_nocase$COHRAID

ECC16S_month24 <- subset_samples(ECC_16S, COHRAID %in% samples_visit24_nocase$COHRAID)
sample_data(ECC16S_month24) <- samples_visit24_nocase

saveRDS(ECC16S_month24, file=file.path(folder, "16S/phyloseq/phyasv_visit24.rds"))


