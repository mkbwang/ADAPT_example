# clean data and retain counts on the phylum level

rm(list=ls())

library(readxl)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(dirmult)

input_folder <- 'real_data/global_gut/'
output_folder <- "simulation/"


# metadata read in 2

meta_data=read_excel(file.path(input_folder, "global_gut_metadata.xls"), sheet = 1, skip = 2) %>%
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
otu_table <- read_tsv(file.path(input_folder, "global_gut_count.txt")) |> as.data.frame()
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
# global_gut_adult <- list(otu_table=otu_table_adult,
#                          taxonomy_table = taxonomies,
#                          metadata=meta_data_adult)
# 
# saveRDS(global_gut_adult, file="real_data/global_gut/global_gut_adult.rds")



# subset of US adults
meta_data_USadult <- meta_data_adult %>% filter(country=="USA")
num_USadult <- nrow(meta_data_USadult) 
otu_table_USadult <- otu_table_adult[, rownames(meta_data_USadult)]
prevalence_USadult <- rowMeans(otu_table_USadult != 0)
set.seed(2023)
qualified_taxa <- which(prevalence_USadult > 0.2)
subset_taxa <- sample(qualified_taxa, size=1000)

otu_table_USadult <- otu_table_USadult[subset_taxa, ]
taxonomies_USadult <- taxonomies[subset_taxa, ]
USadult_metadata <- data.frame(X=c(rep(0, (num_USadult-1)/2), rep(1, (num_USadult+1)/2)))

template <- list(otu_table=otu_table_USadult,
                 metadata=USadult_metadata,
                 taxonomy=taxonomies_USadult)

saveRDS(template, file.path(output_folder, 'global_gut_USadult.rds'))


dirmult_fit <- dirmult(data=t(otu_table_USadult))
gamma_estimate <- sort(dirmult_fit$gamma, decreasing=T)

saveRDS(gamma_estimate, file.path(output_folder, "DM_param.rds"))

