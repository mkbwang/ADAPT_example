library(SQMtools)
library(dplyr)
library(phyloseq)

folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC"
# metadata
load(file.path(folder, "metadata", "metag_metashort.Rdata"))


# Plaque SQM
plaque_object = loadSQM(file.path(folder, "WGS", "metasqueeze", "plaque_assembly"), 
                        engine = 'data.table')
species_counts <- plaque_object$taxa$species$abund
taxonomies <- plaque_object$misc$tax_names_long$species
taxonomy_table <- sapply(taxonomies, 
                         function(longname) strsplit(longname, split=";")[[1]] )
taxonomy_table <- t(taxonomy_table)
taxon_filter <- taxonomy_table[, 1] == "k_Bacteria"
filtered_taxonomy_table <- taxonomy_table[taxon_filter, ]
colnames(filtered_taxonomy_table) <- c("Kingdom", "Phylum", "Class",
                                       "Order", "Family", "Genus", "Species")
species_counts_subset <- species_counts[taxon_filter, ]


plaque_metadata <- metag_metashort %>% filter(SampleType == "plaque")
rownames(plaque_metadata) <- plaque_metadata$OriginalID
plaque_metadata$iscase <- plaque_metadata$CaseStatus == "Case"
phy_metag_plaque <- phyloseq(otu_table(species_counts_subset, taxa_are_rows = T),
                             tax_table(taxonomy_table),
                             sample_data(plaque_metadata))


saveRDS(phy_metag_plaque, file=file.path(folder, "WGS", "phy_metag_plaque.rds"))


