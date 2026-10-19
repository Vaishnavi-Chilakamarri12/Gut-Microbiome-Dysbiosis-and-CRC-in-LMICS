
BiocManager::install("phyloseq")
install.packages(c("ggplot2", "FSA"))
library(phyloseq)
library(ggplot2)
library(FSA)

load("D:/R workflow/DADA2_Rfiles/Kenya_DADA2.rData")
load("D:/R workflow/DADA2_Rfiles/Sri_Lanka_DADA2.rData")
load("D:/R workflow/DADA2_Rfiles/USA_DADA2.rData")
load("D:/R workflow/microbial_diversity.rData")

# Part 1. Calculate alpha diversity for all samples
alpha_all = estimate_richness(physeq, measures = c("Shannon", "Simpson", "Observed"))

# Add Sample ID
alpha_all$SampleID = rownames(alpha_all)

# Add metadata
alpha_all = merge(alpha_all, as.data.frame(sample_data(physeq)),
  by.x = "SampleID",by.y = "row.names")

# Part 2. Define population groups
alpha_all$Population = with(
  alpha_all,
  ifelse(
    Country == "Sri Lanka",
    "Sri Lanka",
    ifelse(
      Country == "Kenya",
      "Kenya",
      ifelse(
        Population_background == "African American",
        "African American",
        ifelse(
          Population_background == "White American",
          "White American",
          NA
        )
      )
    )
  )
)

# Check population assignments
table(alpha_all$Population, useNA = "ifany")
table(alpha_all$Population, alpha_all$Condition) # Check sample counts by population and condition

nsamples(physeq)
sample_names(physeq)
table(sample_data(physeq)$Country)
table(
  sample_data(physeq)$Country,
  sample_data(physeq)$Condition
)

# part 3. SRI LANKA: CRC vs independent healthy controls
# UNPAIRED WILCOXON RANK-SUM TEST
alpha_SL = subset(
  alpha_all,
  Population == "Sri Lanka"
)

# Shannon
wilcox.test(
  Shannon ~ Condition,
  data = alpha_SL,
  exact = FALSE
)

# Simpson
wilcox.test(
  Simpson ~ Condition,
  data = alpha_SL,
  exact = FALSE
)

# Observed richness
wilcox.test(
  Observed ~ Condition,
  data = alpha_SL,
  exact = FALSE
)

# part 4. KENYA: CRC vs independent healthy controls
# UNPAIRED WILCOXON RANK-SUM TEST
alpha_Kenya = subset(
  alpha_all,
  Population == "Kenya"
)

# Shannon
wilcox.test(
  Shannon ~ Condition,
  data = alpha_Kenya,
  exact = FALSE
)

# Simpson
wilcox.test(
  Simpson ~ Condition,
  data = alpha_Kenya,
  exact = FALSE
)

# Observed richness
wilcox.test(
  Observed ~ Condition,
  data = alpha_Kenya,
  exact = FALSE
)


# part 5: USA AFRICAN-AMERICAN COHORT: CRC tissue vs matched normal mucosa
# Load USA pairing metadata
USA_metadata = read.csv("D:/summer research MSU/SRA_data/USA_Metadata.csv", check.names = TRUE, stringsAsFactors = FALSE)

# PAIRED WILCOXON SIGNED-RANK TEST
# Add USA pairing information to alpha_all
USA_pairing = USA_metadata[, c("SRA.number", "Sample_ID")]

colnames(USA_pairing)[1] = "SampleID"

alpha_all = merge(
  alpha_all,
  USA_pairing,
  by = "SampleID",
  all.x = TRUE
)

# AFRICAN AMERICAN
alpha_AA = subset(
  alpha_all,
  Population == "African American"
)

table(alpha_AA$Condition) # Check sample counts
table(alpha_AA$Sample_ID, alpha_AA$Condition) # Check the pair identifier


# Paired Wilcoxon
# Shannon
AA_Shannon = reshape(
  alpha_AA[, c("Sample_ID", "Condition", "Shannon")],
  idvar = "Sample_ID",
  timevar = "Condition",
  direction = "wide"
)

wilcox.test(
  AA_Shannon$Shannon.CRC,
  AA_Shannon$Shannon.Control,
  paired = TRUE,
  exact = FALSE
)


# Simpson
AA_Simpson = reshape(
  alpha_AA[, c("Sample_ID", "Condition", "Simpson")],
  idvar = "Sample_ID",
  timevar = "Condition",
  direction = "wide"
)

wilcox.test(
  AA_Simpson$Simpson.CRC,
  AA_Simpson$Simpson.Control,
  paired = TRUE,
  exact = FALSE
)


# Observed richness
AA_Observed = reshape(
  alpha_AA[, c("Sample_ID", "Condition", "Observed")],
  idvar = "Sample_ID",
  timevar = "Condition",
  direction = "wide"
)

wilcox.test(
  AA_Observed$Observed.CRC,
  AA_Observed$Observed.Control,
  paired = TRUE,
  exact = FALSE
)

# part 6. U.S. WHITE AMERICAN COHORT: CRC tissue vs matched normal mucosa
alpha_WA = subset(
  alpha_all,
  Population == "White American"
)

table(alpha_WA$Condition) # Check sample counts
table(alpha_WA$Sample_ID, alpha_WA$Condition) # Check that each Sample_ID has one CRC and one Control

# Shannon
WA_Shannon = reshape(
  alpha_WA[, c("Sample_ID", "Condition", "Shannon")],
  idvar = "Sample_ID",
  timevar = "Condition",
  direction = "wide"
)

wilcox.test(
  WA_Shannon$Shannon.CRC,
  WA_Shannon$Shannon.Control,
  paired = TRUE,
  exact = FALSE
)


# Simpson
WA_Simpson = reshape(
  alpha_WA[, c("Sample_ID", "Condition", "Simpson")],
  idvar = "Sample_ID",
  timevar = "Condition",
  direction = "wide"
)

wilcox.test(
  WA_Simpson$Simpson.CRC,
  WA_Simpson$Simpson.Control,
  paired = TRUE,
  exact = FALSE
)


# Observed richness
WA_Observed = reshape(
  alpha_WA[, c("Sample_ID", "Condition", "Observed")],
  idvar = "Sample_ID",
  timevar = "Condition",
  direction = "wide"
)

wilcox.test(
  WA_Observed$Observed.CRC,
  WA_Observed$Observed.Control,
  paired = TRUE,
  exact = FALSE
)

####  KRUSKAL-WALLIS TEST  ####
# Compare populations within CRC only
alpha_crc = subset(
  alpha_all,
  Condition == "CRC" &
    !is.na(Population)
)

# Shannon
kruskal.test(
  Shannon ~ Population,
  data = alpha_crc
)

# Simpson
kruskal.test(
  Simpson ~ Population,
  data = alpha_crc
)

# Observed richness
kruskal.test(
  Observed ~ Population,
  data = alpha_crc
)

####  DUNN'S POST-HOC TEST  ####
 
# Shannon
dunnTest(
  Shannon ~ Population,
  data = alpha_crc,
  method = "bh"
)

# Simpson
dunnTest(
  Simpson ~ Population,
  data = alpha_crc,
  method = "bh"
)

# Observed richness
dunnTest(
  Observed ~ Population,
  data = alpha_crc,
  method = "bh"
)

#save.image(file = "Alphadiversity_statistical.rData")
