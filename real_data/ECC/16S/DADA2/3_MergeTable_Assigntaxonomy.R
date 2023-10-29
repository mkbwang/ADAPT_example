library(dada2); packageVersion("dada2")


folder <- "/nfs/turbo/sph-ligen/wangmk/ADAPT_example/real_data/ECC/16S"
# Merge multiple runs 
table1 <- readRDS(file.path(folder, 'seqfetch', 'filter_sequences', 'seqtable_1.rds'))
table2 <- readRDS(file.path(folder, 'seqfetch', 'filter_sequences', 'seqtable_2.rds'))
table3 <- readRDS(file.path(folder, 'seqfetch', 'filter_sequences', 'seqtable_3.rds'))
table_all <- mergeSequenceTables(table1, table2, table3)

# Remove chimeras
seqtab <- removeBimeraDenovo(table_all, method="consensus", multithread=TRUE)

# Assign taxonomy
## reference 16s rRNA genomes
taxa_toSp <- assignTaxonomy(seqtab, 
                            file.path(folder, 'DADA2', "HOMD_AssignTaxaToSpecies_plusextras.fasta"), multithread = TRUE)

taxa_subset <- taxa_toSp[!is.na(taxa_toSp[, 7]), ]


# Write to disk
saveRDS(seqtab, file.path(folder, 'phyloseq', "seqtab_homd.rds")) 
saveRDS(taxa_toSp, file.path(folder, 'phyloseq', "tax_toSp_homd.rds")) 

