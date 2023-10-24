# clean data and retain counts on the phylum level

rm(list=ls())

library(readxl)
library(tidyverse)
library(dplyr)
library(nloptr)
library(phyloseq)
library(stringr)
library(ggpubr)
library(magrittr)

folder <- '/scratch/ligen_root/ligen0/wangmk/POLDA_example'
output_folder <- "/scratch/ligen_root/ligen0/wangmk/POLDA_example"


# metadata read in 2

meta_data=read_excel(file.path(folder, "real_data/global_gut/NIHMS365354-supplement-3.xls"), sheet = 1, skip = 2) %>%
  as.data.frame()
meta_data=meta_data[!is.na(meta_data$`Sample Identifier`), ]
meta_data=meta_data%>%
  transmute(Sample.ID=`Sample Identifier`, age=`Age (Years)`, gender=Gender,country=Country,
            Depth = `Number of V4-16S rRNA Illumina sequences`)%>%
  arrange(Sample.ID)
meta_data$age=signif(as.numeric(meta_data$age), digits = 2)



# read in taxonomy
# taxonomy=read_tsv("real_data/Yatsunenko_data/global_gut_taxonomy.txt") |> as.data.frame()
# rownames(taxonomy) <- as.character(taxonomy$OTU_ID)
# taxonomy$OTU_ID <- NULL


# read in the OTU table and aggregate into phylum level
otu_table <- read_tsv(file.path(folder, "real_data/global_gut/Yatsunenko2012.txt")) |> as.data.frame()
rownames(otu_table) <- as.character(otu_table$OTU_ID)


taxonomy <- otu_table[, 530]
taxonomies <- sapply(taxonomy, function(fullname) strsplit(fullname, split=";")[[1]]) |>
  t()
rownames(taxonomies) <- rownames(otu_table)
colnames(taxonomies) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# taxonomies[, 6] <- str_replace(taxonomies[, 6], pattern = "g__", replacement = "")

otu_table <- otu_table[, -c(1, 530)] # get rid of taxonomy column

# otu_table$Genus <- taxonomies[taxonomies[, 6] != "", 6]
# otu_table <- otu_table %>% group_by(Genus) %>% summarize(across(everything(), sum))
# otu_table <- otu_table[, -1]

## the subset of adults
meta_data_adult <- meta_data %>% filter(age >= 18) %>% as.data.frame()
rownames(meta_data_adult) <- meta_data_adult$Sample.ID
meta_data_adult$Sample.ID <- NULL
otu_table_adult <- otu_table[, rownames(meta_data_adult)] |> as.matrix()
global_gut_adult <- list(otu_table=otu_table_adult,
                         taxonomy_table = taxonomies,
                         metadata=meta_data_adult)

saveRDS(global_gut_adult, file="real_data/global_gut/global_gut_adult.rds")



# subset of US adults
meta_data_USadult <- meta_data_adult %>% filter(country=="USA")
num_USadult <- nrow(meta_data_USadult) 
otu_table_USadult <- otu_table_adult[, rownames(meta_data_USadult)]
prevalence_USadult <- rowMeans(otu_table_USadult != 0)
set.seed(2023)
qualified_taxa <- which(prevalence_USadult > 0.1)
subset_taxa <- sample(qualified_taxa, size=1000)

otu_table_USadult <- otu_table_USadult[subset_taxa, ]
taxonomies_USadult <- taxonomies[subset_taxa, ]
USadult_metadata <- data.frame(X=c(rep(0, (num_USadult-1)/2), rep(1, (num_USadult+1)/2)))

template <- list(otu_table=otu_table_USadult,
                 metadata=USadult_metadata,
                 taxonomy=taxonomies_USadult)

saveRDS(template, file.path(output_folder, 'real_data/global_gut/global_gut_USadult.rds'))

library(dirmult)

dirmult_fit <- dirmult(data=t(otu_table_USadult))

gamma_estimate <- sort(dirmult_fit$gamma)

saveRDS(gamma_estimate, file.path(output_folder, "simulation/simulate_DirMultinom/global_gut_template_dirparam.rds"))

# ## the subset of infants
# meta_data_infant <- meta_data %>% filter(age <=0.5)
# rownames(meta_data_infant) <- meta_data_infant$Sample.ID
# meta_data_infant$Sample.ID <- NULL
# otu_table_infant <- otu_table[, rownames(meta_data_infant)] |> as.matrix()
# # prevalence_infant <- rowMeans(otu_table_infant != 0)
# # 
# # otu_table_infant <- otu_table_infant[prevalence_infant > 0.2, ]
# # taxonomies_infant <- taxonomies[prevalence_infant > 0.2, ]
# global_gut_infant <- list(otu_table=otu_table_infant,
#                           taxonomy_table = taxonomies,
#                           metadata=meta_data_infant)
# 
# saveRDS(global_gut_infant, file="real_data/global_gut/global_gut_infant.rds")
