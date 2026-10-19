# Install packages
BiocManager::install("phyloseq")
install.packages("vegan") # For Bray-Curtis
install.packages("ggplot2")
library(phyloseq)
library(vegan)
library(ggplot2)

# Load the data
load("D:/R workflow/DADA2_Rfiles/DADA2_kenya.rData")
load("D:/R workflow/DADA2_Rfiles/DADA2_Srilanka.rData")
load("D:/R workflow/DADA2_Rfiles/DADA2_USA.rData")
physeq = readRDS("D:/R workflow/physeq_all.rData")

metadata = read.csv("D:/summer research MSU/SRA_data/DADA2_Metadata.csv", row.names = 1, check.names = TRUE)
table(metadata$Population_background)
sample_data(physeq)$Condition = metadata[rownames(sample_data(physeq)), "Condition"]

# Creating two separate phyloseq objects: one for CRC patients and one for controls.
physeq_crc = subset_samples(physeq, Condition == "CRC")
physeq_n  = subset_samples(physeq, Condition == "Control")

# Number of taxa that occur exactly once across all samples:
#sum(taxa_sums(physeq) == 1)
#table(taxa_sums(physeq))  # distribution of taxa abundances

# CRC: Estimate alpha diversity metrics
alpha_crc = estimate_richness(physeq_crc, measures = c("Shannon", "Simpson", "Observed"))
# Add metadata back for grouping
alpha_crc$SampleID = rownames(alpha_crc)
alpha_crc = merge(alpha_crc, as.data.frame(sample_data(physeq_crc)), 
                  by.x = "SampleID", by.y = "row.names")
# Plot Shannon diversity
ggplot(alpha_crc, aes(x = Group, y = Shannon, fill = Country)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Shannon Diversity Index - CRC", y = "Shannon Index", x = "Group")

# Control: Estimate alpha diversity metrics
alpha_n = estimate_richness(physeq_n, measures = c("Shannon", "Simpson", "Observed"))
# Add metadata back for grouping
alpha_n$SampleID = rownames(alpha_n)
alpha_n = merge(alpha_n, as.data.frame(sample_data(physeq_n)), 
                  by.x = "SampleID", by.y = "row.names")
# Plot Shannon diversity for Control samples
ggplot(alpha_n, aes(x = Group, y = Shannon, fill = Country)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Shannon Diversity Index - Control", y = "Shannon Index", x = "Group")

################ ALPHA DIVERSITY: EACH COUNTRY ################
## Notes: This workflow is reused for each country or Population_background by changing only the subsetting variable (e.g., Country or Population_background) and the group labels

# Subset African American samples only
physeq_usa = subset_samples(
  physeq,
  Population_background == "African American"
)

# Estimate alpha diversity metrics (Shannon and Simpson)
alpha_usa = estimate_richness(
  physeq_usa,
  measures = c("Shannon", "Simpson")
)

# Merge metadata for plotting
alpha_usa$SampleID = rownames(alpha_usa)

alpha_usa = merge(
  alpha_usa,
  as.data.frame(sample_data(physeq_usa)),
  by.x = "SampleID",
  by.y = "row.names"
)

# Make Group a factor with CRC first, Control second
alpha_usa$Group = factor(
  alpha_usa$Group,
  levels = c("AA_CRC", "AA_N"),
  labels = c("CRC", "Control")
)

# Plot Shannon diversity
ggplot(alpha_usa, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "USA: African American",
    y = "Shannon Index",
    x = "Condition"
  ) +
  scale_fill_manual(
    values = c("CRC" = "tomato", "Control" = "skyblue")
  )

############### BETA DIVERSITY: CRC ONLY ###############

# Subset CRC samples
physeq_crc = subset_samples(physeq, Condition == "CRC")
physeq_crc_rel = transform_sample_counts(physeq_crc, function(x) x / sum(x)) #Transform to relative abundance

# Bray-Curtis distance matrix
bray_crc = distance(physeq_crc_rel, method = "bray")
ord_crc = ordinate(physeq_crc_rel, method = "PCoA", distance = bray_crc) #PCoA ordination

# Plot pcoA for CRC
plot_ordination(physeq_crc_rel, ord_crc, color = "Population_background") +
  geom_point(size = 4, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Bray-Curtis Distance: CRC samples", color = "Population_background")

#save.image(file = "microbial_diversity.rData")
