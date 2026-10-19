# Gut Microbiome Dysbiosis and CRC in Low-and-Middle-Income countries

## Overview
This repository contains Metadata and R scirpts that were used in this cross-population study of colorectal cancer (CRC)-associated mucosal microbiomes. The study used publicly available microbiome data from Sri Lanka, Kenya, and the United States (USA). Standardized bioinformatics pipelines were used to perform taxonomic profiling, microbial diversity analysis and statistical testing.

## Study Workflow
### 1. Data Collection
- Systematic literature review.
- Retrieval of publicly available sequencing data from NCBI SRA and EMBL-EBI.
- sample metadata compiled from public database.
### 2. Data Processing
- Downloading raw paired-end FASTQ files using SRA Toolkit.
- DADA2 analysis (Quality control, filtering, trimming, denoising, and chimera removal).
- ASV inference and taxonomic classification using the SILVA reference database (release 138.2).
- Integration of ASV tables, taxonomy, and metadata into phyloseq objects.
### 3. Microbiome Analysis
- Taxonomic composition analysis at the phylum level.
- Alpha diversity analysis using Shannon, Simpson, and observed ASV richness.
- Beta diversity analysis using Bray–Curtis dissimilarity and PCoA.
### 4. Statistical Analysis
- Comparison of CRC and control groups using Wilcoxon rank-sum and signed-rank tests, as appropriate.
- Comparison of CRC samples across populations using Kruskal–Wallis tests and Dunn’s post-hoc comparisons.
- Evaluation of microbial community differences using PERMANOVA.
