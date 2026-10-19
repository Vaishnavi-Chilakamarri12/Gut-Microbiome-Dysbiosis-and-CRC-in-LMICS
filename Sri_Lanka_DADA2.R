# Install necessary packages
install.packages("BiocManager")
BiocManager::install("dada2", ask = FALSE)
install.packages("Rcpp")
library(Rcpp)
library(dada2)
library(ShortRead)
library(stats)

#Define paths for SL data
path_CRC_SL = "D:/summer research MSU/SRA_data/Srilanka/SL_CRC"
path_N_SL = "D:/summer research MSU/SRA_data/Srilanka/SL_N"
#list.files(path_CRC_SL)
#list.files(path_N_SL)

#Sort and list fastq files
fnFs_CRC = sort(list.files(path_CRC_SL, pattern="_1.fastq", full.names = TRUE))
fnRs_CRC = sort(list.files(path_CRC_SL, pattern="_2.fastq", full.names = TRUE))
fnFs_N  = sort(list.files(path_N_SL, pattern="_1.fastq", full.names = TRUE))
fnRs_N  = sort(list.files(path_N_SL, pattern="_2.fastq", full.names = TRUE))

#Create a Vector to extract sample names up to the first underscore
SLCRC_sample.names = sapply(strsplit(basename(fnFs_CRC), "_"), `[`, 1)
SLN_sample.names = sapply(strsplit(basename(fnFs_N), "_"), `[`, 1)
fnRs_CRC[1:7]
fnFs_N[1:7]

#Inspect read quality profiles
#plotQualityProfile(fnFs_CRC[1:7])
#plotQualityProfile(fnRs_CRC[1:7])
#plotQualityProfile(fnFs_N[1:7])
#plotQualityProfile(fnRs_N[1:7])

#Part 2: Filter and trim data
SL_CRCfiltFs = file.path(path_CRC_SL, "filtered", paste0(SLCRC_sample.names, "_SLCRC_F_filt.fastq.gz"))
SL_CRCfiltRs = file.path(path_CRC_SL, "filtered", paste0(SLCRC_sample.names, "_SLCRC_R_filt.fastq.gz"))
SL_NfiltFs = file.path(path_N_SL, "filtered", paste0(SLN_sample.names, "_SLN_F_filt.fastq.gz"))
SL_NfiltRs = file.path(path_N_SL, "filtered", paste0(SLN_sample.names, "_SLN_R_filt.fastq.gz"))

names(SL_CRCfiltFs) = SLCRC_sample.names
names(SL_CRCfiltRs) = SLCRC_sample.names
names(SL_NfiltFs) = SLN_sample.names
names(SL_NfiltRs) = SLN_sample.names

SLCRC_out = filterAndTrim(fnFs_CRC, SL_CRCfiltFs, fnRs_CRC, SL_CRCfiltRs, 
                            truncLen = c(220, 150), maxN = 0, maxEE=c(2,2), truncQ=2, 
                            rm.phix = TRUE, compress = TRUE, multithread = FALSE)
SLN_out = filterAndTrim(fnFs_N, SL_NfiltFs, fnRs_N, SL_NfiltRs, 
                          truncLen = c(220, 150), maxN = 0, maxEE=c(2,2), truncQ=2, 
                          rm.phix = TRUE, compress = TRUE, multithread = FALSE)
#SLCRC_out
#SLN_out

#part 3: learn error rates
errF_SLCRC = learnErrors(SL_CRCfiltFs, multithread=TRUE)
errR_SLCRC = learnErrors(SL_CRCfiltRs, multithread=TRUE)
errF_SLN = learnErrors(SL_NfiltFs, multithread=TRUE)
errR_SLN = learnErrors(SL_NfiltRs, multithread=TRUE)

# Plot error rates for SL samples
plotErrors(errF_SLCRC, nominalQ=TRUE)
plotErrors(errR_SLCRC, nominalQ = TRUE)
plotErrors(errF_SLN, nominalQ=TRUE)
plotErrors(errR_SLN, nominalQ = TRUE)
#black line-estimated error rates after convergence of the machine-learning algorithm
#red line- error rates expected under the nominal definition of the Q-score
#Here the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected.
##estimate error rates of African American associated with sequencing process

#part 4: ASV Sample inference (used to filter the trimmed sequence data)
SL_CRCdadaFs = dada(SL_CRCfiltFs, err=errF_SLCRC, multithread=TRUE)
SL_CRCdadaRs = dada(SL_CRCfiltRs, err=errR_SLCRC, multithread=TRUE)
SL_NdadaFs = dada(SL_NfiltFs, err=errF_SLN, multithread=TRUE)
SL_NdadaRs = dada(SL_NfiltRs, err=errR_SLN, multithread=TRUE)

#part 5: Merge forward and reverse reads
SL_CRCmergers = mergePairs(SL_CRCdadaFs, SL_CRCfiltFs, SL_CRCdadaRs, SL_CRCfiltRs, verbose=TRUE)
SL_Nmergers = mergePairs(SL_NdadaFs, SL_NfiltFs, SL_NdadaRs, SL_NfiltRs, verbose=TRUE)

#part 6: Construction of amplicon sequence variant (ASV) table
SL_CRCseqtab = makeSequenceTable(SL_CRCmergers)
SL_Nseqtab = makeSequenceTable(SL_Nmergers)

# Part 7: removing chimeras
SL_CRCseqtab.nochim = removeBimeraDenovo(SL_CRCseqtab, method="consensus",verbose=TRUE)
table(nchar(getSequences(SL_CRCseqtab.nochim)))
sum(SL_CRCseqtab.nochim)/sum(SL_CRCseqtab)
SL_Nseqtab.nochim = removeBimeraDenovo(SL_Nseqtab, method="consensus", verbose=TRUE)
table(nchar(getSequences(SL_Nseqtab.nochim)))
sum(SL_Nseqtab.nochim)/sum(SL_Nseqtab)

rownames(SL_CRCseqtab.nochim)
rownames(SL_Nseqtab.nochim)

#part 8: Tracking the reads
getN = function(x) sum(getUniques(x))
SL_CRCtrack = cbind(SLCRC_out, sapply(SL_CRCdadaFs, getN), sapply(SL_CRCdadaRs, getN), 
                     sapply(SL_CRCmergers, getN), rowSums(SL_CRCseqtab.nochim))
SL_Ntrack = cbind(SLN_out, sapply(SL_NdadaFs, getN), sapply(SL_NdadaRs, getN), 
                    sapply(SL_Nmergers, getN), rowSums(SL_Nseqtab.nochim))
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(SL_CRCtrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(SL_CRCtrack) = SLCRC_sample.names
colnames(SL_Ntrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(SL_Ntrack) = SLN_sample.names

#part 9: assigning taxonomy using SILVA database
SL_CRCtaxa = assignTaxonomy(SL_CRCseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
SL_Ntaxa = assignTaxonomy(SL_Nseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#inspect the taxonomic assignments
SL_CRCtaxa.print = SL_CRCtaxa
rownames(SL_CRCtaxa.print) = NULL
head(SL_CRCtaxa.print)
colnames(SL_CRCtaxa)


#save.image(file = "DADA2_Srilanka.rData")