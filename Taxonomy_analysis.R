# Install and load required packages
BiocManager::install(c("phyloseq", "dada2"), force = TRUE)
install.packages("Rcpp")
install.packages("ggplot2")
library(phyloseq)
library(ggplot2)
library(Rcpp)
library(dada2)
library(dplyr)
library(tidyr)

# Load the data, tools and metadata
load("D:/R workflow/DADA2_Rfiles/Kenya_DADA2.rData")
load("D:/R workflow/DADA2_Rfiles/Sri_Lanka_DADA2.rData")
load("D:/R workflow/DADA2_Rfiles/USA_DADA2.rData")

silva = "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz"
gtdb_species = "D:/summer research MSU/Tools/GTDB_dada2_assignment_species.fa.gz"
metadata = read.csv("D:/summer research MSU/SRA_data/DADA2_Metadata.csv", check.names = FALSE, row.names = 1)

# Combine CRC and HC ASV tables of all countries
SL_seqtab_combined = mergeSequenceTables(SL_CRCseqtab.nochim, SL_Nseqtab.nochim)
ken_seqtab_combined = mergeSequenceTables(ken_CRCseqtab.nochim, ken_Nseqtab.nochim)
USA_seqtab_combined = mergeSequenceTables(WA_CRCseqtab.nochim, WA_Nseqtab.nochim, 
                                          AA_CRCseqtab.nochim, AA_Nseqtab.nochim)

# Merge all ASV, taxonomy tables and add species table
seqtab_combined = mergeSequenceTables(SL_seqtab_combined, USA_seqtab_combined, ken_seqtab_combined)
tax_combined = assignTaxonomy(seqtab_combined, silva, multithread = TRUE)
tax_combined_species = addSpecies(tax_combined, gtdb_species)

# Create phyloseq objects
otu = otu_table(seqtab_combined, taxa_are_rows = FALSE)
tax = tax_table(tax_combined_species)
metadata = metadata[rownames(otu), ] #to ensure sample names match
sam = sample_data(metadata)

sample_names(otu)
rownames(metadata)

# Clean metadata
clean_names = gsub("\t", "", rownames(metadata)) #remove tab characters
clean_names = trimws(clean_names) #Trim white spaces
rownames(metadata) = clean_names #Reassign cleaned names

# Final phyloseq object
physeq = phyloseq(otu, tax, sam)

#head(tax_table(physeq))
#saveRDS(physeq, "D:/R workflow/physeq_all.rData")

############  TAXONOMIC BARPLOT  ###################
# Agglomerate at Phylum level
physeq_phylum = tax_glom(physeq, taxrank = "Phylum")

# Subset by condition
physeq_crc = subset_samples(physeq_phylum, Condition == "CRC")
physeq_n  = subset_samples(physeq_phylum, Condition == "Control")

# Merge samples by country for each condition
physeq_crc_merged = merge_samples(physeq_crc, group = "Group")
physeq_n_merged  = merge_samples(physeq_n, group = "Group")

# Transform to relative abundance
physeq_crc_rel = transform_sample_counts(physeq_crc_merged, function(x) x / sum(x))
physeq_n_rel  = transform_sample_counts(physeq_n_merged, function(x) x / sum(x))

# Plot CRC barplot
plot_bar(physeq_crc_rel, fill = "Phylum") +
  ggtitle("Taxonomic composition of CRC samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Group") +
  labs(fill = "Phylum")

# Plot control barplot
plot_bar(physeq_n_rel, fill = "Phylum") +
  ggtitle("Taxonomic composition of control samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Group") +
  labs(fill = "Phylum")

save.image(file = "taxonomic_barplot.rData")

########### Percentage (%) of major phyla by cohort #############
# remove taxa (phyla) that have zero total counts across all samples
physeq_phylum = prune_taxa(taxa_sums(physeq_phylum) > 0, physeq_phylum)

#Transform to relative abundance
physeq_rel = transform_sample_counts(physeq_phylum, function(x) x / sum(x))

#Convert to long dataframe
phylum_df = psmelt(physeq_rel)  # Sample, Group, Condition, Phylum, Abundance

#Calculating mean % relative abundance per Phylum per Country and Condition
phylum_summary = phylum_df %>%
  group_by(Group, Condition, Phylum) %>%
  summarise(Mean_RA = mean(Abundance) * 100, SD_RA   = sd(Abundance) * 100, .groups = "drop")

#Creating combined label for columns
phylum_summary = phylum_summary %>%
  mutate(Group_Condition = paste(Group, Condition, sep = "_")) %>%
  select(Group_Condition, Phylum, Mean_RA,SD_RA)

#Pivot wider: Phylum as rows, Group_Condition as columns
phylum_table = phylum_summary %>%
  pivot_wider(names_from = Group_Condition, values_from = c(Mean_RA, SD_RA))

#write.csv(phylum_table, "D:/R workflow/Phylum_Percent.csv", row.names = TRUE)
