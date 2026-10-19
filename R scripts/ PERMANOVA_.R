#######################  PERMANOVA  ############################
BiocManager::install("phyloseq")
install.packages("vegan")
install.packages("remotes")
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("cluster")
install.packages("permute")

library(permute)
library(cluster)
library(pairwiseAdonis)
library(phyloseq)
library(vegan)

physeq = readRDS("D:/R workflow/physeq_all.rData") #load data

#step 1. Subset samples (CRC + Control only)
physeq_all = subset_samples(
  physeq,
  Condition %in% c("CRC", "Control")
)

# To Prevent issues in distance calculations, remove taxa with zero counts
physeq_all = prune_taxa(taxa_sums(physeq_all) > 0, physeq_all)

# step 2. Bray–Curtis distance calculation
bray_all = phyloseq::distance(physeq_all, method = "bray")

#step 3. Metadata preparation
meta_all = data.frame(sample_data(physeq_all))
meta_all = meta_all[labels(bray_all), ]  # This step is to ensure correct ordering
#Convert variables to factors
meta_all$Condition = factor(meta_all$Condition)
meta_all$Country   = factor(meta_all$Country)
meta_all$Group     = factor(meta_all$Group)
meta_all$Population_background = factor(meta_all$Population_background)

#step 4. PERMANOVA Part A: CRC vs control within each Country
cat("\n============================\n")
cat("PERMANOVA: CRC vs control within each Country\n")
cat("============================\n")

countries = unique(meta_all$Country)
perm_results_within = list()

#Loop through each country
for (ctry in countries) {           
  
  cat("\n--- Country:", ctry, "---\n")
  
  meta_sub = subset(meta_all, Country == ctry)
  
  bray_sub = as.dist(
    as.matrix(bray_all)[rownames(meta_sub), rownames(meta_sub)]
  )
  
  # PERMANOVA: Condition effect
  adonis_res = adonis2(
    bray_sub ~ Condition,
    data = meta_sub,
    permutations = 999
  )
  
  print(adonis_res)
  
  # PERMDISP
  cat("\nDispersion test (PERMDISP):\n")
  bd = betadisper(bray_sub, meta_sub$Condition)
  disp_res = permutest(bd, permutations = 999)
  
  print(disp_res)
  
  perm_results_within[[as.character(ctry)]] = list(
    permanova = adonis_res,
    dispersion = disp_res
  )
}

#step 5. PERMANOVA Part B: Global test across all Group strata
cat("PERMANOVA: All cohort-condition groups (Group)\n")

adonis_group = adonis2(
  bray_all ~ Group,
  data = meta_all,
  permutations = 999
)

print(adonis_group)

# PERMDISP
cat("\nDispersion test (PERMDISP) across Groups:\n")
bd_group = betadisper(bray_all, meta_all$Group)
disp_group = permutest(bd_group, permutations = 999)
print(disp_group)

#step 6. PERMANOVA Part C: CRC-only comparison across populations
cat("PERMANOVA: CRC-only comparison across Countries\n")

meta_crc = subset(meta_all, Condition == "CRC")
bray_crc = as.dist(
  as.matrix(bray_all)[rownames(meta_crc), rownames(meta_crc)]
)

adonis_crc_country = adonis2(
  bray_crc ~ Country,
  data = meta_crc,
  permutations = 999
)
print(adonis_crc_country)

# PERMDISP
cat("\nDispersion test (PERMDISP) across Countries in CRC-only:\n")
bd_crc = betadisper(bray_crc, meta_crc$Country)
disp_crc = permutest(bd_crc, permutations = 999)
print(disp_crc)

#step 7. USA CRC-only Population_background comparison (AA vs WA)

cat("PERMANOVA: Population_background within USA CRC-only\n")

meta_usa_crc = subset(meta_crc, Country == "USA")

if (length(unique(meta_usa_crc$Population_background)) > 1) {
  
  bray_usa_crc = as.dist(
    as.matrix(bray_all)[rownames(meta_usa_crc), rownames(meta_usa_crc)]
  )
  
  adonis_usa_eth = adonis2(
    bray_usa_crc ~ Population_background,
    data = meta_usa_crc,
    permutations = 999
  )
  
  print(adonis_usa_eth)
  
  # PERMDISP
  cat("\nDispersion test (PERMDISP) within USA CRC-only:\n")
  bd_usa = betadisper(bray_usa_crc, meta_usa_crc$Population_background)
  disp_usa = permutest(bd_usa, permutations = 999)
  print(disp_usa)
  
} else {
  cat("Not enough Population_background groups in USA subset.\n")
}

save.image(file = "PERMANOVA.rData")

#### NOTES ####
#adonis2() in PERMANOVA tests whether groups are different from each other in overall community composition using a distance matrix. It uses permutations (not normal assumptions) and can be affected if groups have very different levels of variation
# betadisper() in PERMDISP Tests whether groups have similar levels of variation (spread) in multivariate space. It should not be significant if you want to confidently interpret PERMANOVA results.
