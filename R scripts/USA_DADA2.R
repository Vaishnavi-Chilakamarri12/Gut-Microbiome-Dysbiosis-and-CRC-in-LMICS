#install packages
install.packages("BiocManager")
BiocManager::install("dada2", force = TRUE)
install.packages("Rcpp")
library(Rcpp)
library(dada2)
library(ShortRead)
library(stats)
#Paths to test and Control
path_CRC_WAUSA = "D:/summer research MSU/SRA_data/USA/White_Americans/WA_CRC"
path_N_WAUSA = "D:/summer research MSU/SRA_data/USA/White_Americans/WA_N"
path_CRC_AAUSA = "D:/summer research MSU/SRA_data/USA/African_americans/AA_CRC"
path_N_AAUSA = "D:/summer research MSU/SRA_data/USA/African_americans/AA_N"
#list.files(path_CRC_AAUSA)
#list.files(path_N_WAUSA)

#Sort and list fastq files
#White Americans
fnFs_CRC_WAUSA = sort(list.files(path_CRC_WAUSA, pattern="_1.fastq", full.names = TRUE))
fnFs_N_WAUSA = sort(list.files(path_N_WAUSA, pattern="_1.fastq", full.names = TRUE))
fnRs_CRC_WAUSA = sort(list.files(path_CRC_WAUSA, pattern="_2.fastq", full.names = TRUE))
fnRs_N_WAUSA = sort(list.files(path_N_WAUSA, pattern="_2.fastq", full.names = TRUE))

#African Americans
fnFs_CRC_AAUSA = sort(list.files(path_CRC_AAUSA, pattern="_1.fastq", full.names = TRUE))
fnFs_N_AAUSA = sort(list.files(path_N_AAUSA, pattern="_1.fastq", full.names = TRUE))
fnRs_CRC_AAUSA = sort(list.files(path_CRC_AAUSA, pattern="_2.fastq", full.names = TRUE))
fnRs_N_AAUSA = sort(list.files(path_N_AAUSA, pattern="_2.fastq", full.names = TRUE))

#Create a Vector to extract sample names up to the first underscore
#White American
WAUSA_sample.names = sapply(strsplit(basename(fnFs_CRC_WAUSA), "_"), `[`, 1)
WAUSAN_sample.names = sapply(strsplit(basename(fnFs_N_WAUSA), "_"), `[`, 1)
WAUSA_sample.names = sapply(strsplit(basename(fnRs_CRC_WAUSA), "_"), `[`, 1)
WAUSAN_sample.names = sapply(strsplit(basename(fnRs_N_WAUSA), "_"), `[`, 1)

#African American
AAUSA_sample.names = sapply(strsplit(basename(fnFs_CRC_AAUSA), "_"), `[`, 1)
AAUSAN_sample.names = sapply(strsplit(basename(fnFs_N_AAUSA), "_"), `[`, 1)
AAUSA_sample.names = sapply(strsplit(basename(fnRs_CRC_AAUSA), "_"), `[`, 1)
AAUSAN_sample.names = sapply(strsplit(basename(fnRs_N_AAUSA), "_"), `[`, 1)

#Inspect read quality profiles
#White American
#plotQualityProfile(fnFs_CRC_WAUSA[1:7])
#plotQualityProfile(fnRs_CRC_WAUSA[1:7])
#plotQualityProfile(fnFs_N_WAUSA[1:7])
#plotQualityProfile(fnRs_N_WAUSA[1:7])

#African American
#plotQualityProfile(fnFs_CRC_AAUSA[1:7])
#plotQualityProfile(fnRs_CRC_AAUSA[1:7])
#plotQualityProfile(fnFs_N_AAUSA[1:7])
#plotQualityProfile(fnRs_N_AAUSA[1:7])

#part 2: filter and trim
#Filter White Americans
WAUSA_CRCfiltFs = file.path(path_CRC_WAUSA, "filtered", paste0(WAUSA_sample.names, "_WAUSACRC_F_filt.fastq.gz"))
WAUSA_CRCfiltRs = file.path(path_CRC_WAUSA, "filtered", paste0(WAUSA_sample.names, "_WACUSACRC_R_filt.fastq.gz"))
WAUSA_NfiltFs = file.path(path_N_WAUSA, "filtered", paste0(WAUSAN_sample.names, "_WAUSAN_F_filt.fastq.gz"))
WAUSA_NfiltRs = file.path(path_N_WAUSA, "filtered", paste0(WAUSAN_sample.names, "_WAUSAN_R_filt.fastq.gz"))
names(WAUSA_CRCfiltFs) = WAUSA_sample.names
names(WAUSA_CRCfiltRs) = WAUSA_sample.names
names(WAUSA_NfiltFs) = WAUSAN_sample.names
names(WAUSA_NfiltRs) = WAUSAN_sample.names
#Trim White Americans
WAUSACRC_out = filterAndTrim(fnFs_CRC_WAUSA, WAUSA_CRCfiltFs, fnRs_CRC_WAUSA, WAUSA_CRCfiltRs, truncLen=c(240,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=FALSE)
WAUSAN_out = filterAndTrim(fnFs_N_WAUSA, WAUSA_NfiltFs, fnRs_N_WAUSA, WAUSA_NfiltRs, truncLen=c(240,200),
                        maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=FALSE)
#WAUSACRC_out
#WAUSAN_out

#Filter African Americans
AAUSA_CRCfiltFs = file.path(path_CRC_AAUSA, "filtered", paste0(AAUSA_sample.names, "_AAUSACRC_F_filt.fastq.gz"))
AAUSA_CRCfiltRs = file.path(path_CRC_AAUSA, "filtered", paste0(AAUSA_sample.names, "_AAUSACRC_R_filt.fastq.gz"))
AAUSA_NfiltFs = file.path(path_N_AAUSA, "filtered", paste0(AAUSAN_sample.names, "_AAUSAN_F_filt.fastq.gz"))
AAUSA_NfiltRs = file.path(path_N_AAUSA, "filtered", paste0(AAUSAN_sample.names, "_AAUSAN_R_filt.fastq.gz"))
names(AAUSA_CRCfiltFs) = AAUSA_sample.names
names(AAUSA_CRCfiltRs) = AAUSA_sample.names
names(AAUSA_NfiltFs) = AAUSAN_sample.names
names(AAUSA_NfiltRs) = AAUSAN_sample.names
#Trim African Americans
AAUSACRC_out = filterAndTrim(fnFs_CRC_AAUSA, AAUSA_CRCfiltFs, fnRs_CRC_AAUSA, AAUSA_CRCfiltRs, truncLen=c(230,180),
                             maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                             compress=TRUE, multithread=FALSE)
AAUSAN_out = filterAndTrim(fnFs_N_AAUSA, AAUSA_NfiltFs, fnRs_N_AAUSA, AAUSA_NfiltRs, truncLen=c(230,180),
                            maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=FALSE)
#AAUSACRC_out
#AAUSAN_out 

#part 3: learn error rates
#estimate error rates of White American associated with sequencing process
errF_WAUSACRC = learnErrors(WAUSA_CRCfiltFs)
errR_WAUSACRC = learnErrors(WAUSA_CRCfiltRs)
errF_WAUSAN = learnErrors(WAUSA_NfiltFs)
errR_WAUSAN = learnErrors(WAUSA_NfiltRs)
#plotErrors(errF_WAUSACRC, nominalQ=TRUE)
#plotErrors(errR_WAUSACRC, nominalQ = TRUE)
#plotErrors(errF_WAUSAN, nominalQ=TRUE)
#plotErrors(errR_WAUSAN, nominalQ = TRUE)

#estimate error rates of African American associated with sequencing process
errF_AAUSACRC = learnErrors(AAUSA_CRCfiltFs)
errR_AAUSACRC = learnErrors(AAUSA_CRCfiltRs)
errF_AAUSAN = learnErrors(AAUSA_NfiltFs)
errR_AAUSAN = learnErrors(AAUSA_NfiltRs)
#plotErrors(errF_AAUSACRC, nominalQ=TRUE)
#plotErrors(errR_AAUSACRC, nominalQ = TRUE)
#plotErrors(errF_AAUSAN, nominalQ=TRUE)
#plotErrors(errR_AAUSAN, nominalQ = TRUE)

#part 4: ASV Sample inference (used to filter the trimmed sequence data)
#White American
WA_CRCdadaFs = dada(WAUSA_CRCfiltFs, err=errF_WAUSACRC)
WA_CRCdadaRs = dada(WAUSA_CRCfiltRs, err=errR_WAUSACRC)
WA_NdadaFs = dada(WAUSA_NfiltFs, err=errF_WAUSAN)
WA_NdadaRs = dada(WAUSA_NfiltRs, err=errR_WAUSAN)
#WA_NdadaFs[[1]] # Access the denoised data for the first sample

#African American
AA_CRCdadaFs = dada(AAUSA_CRCfiltFs, err=errF_AAUSACRC)
AA_CRCdadaRs = dada(AAUSA_CRCfiltRs, err=errR_AAUSACRC)
AA_NdadaFs = dada(AAUSA_NfiltFs, err=errF_AAUSAN)
AA_NdadaRs = dada(AAUSA_NfiltRs, err=errR_AAUSAN)

#part 5: Merge forward and reverse reads
#White American
WA_CRCmergers = mergePairs(WA_CRCdadaFs, WAUSA_CRCfiltFs, WA_CRCdadaRs, WAUSA_CRCfiltRs, verbose=TRUE)
WA_Nmergers = mergePairs(WA_NdadaFs, WAUSA_NfiltFs, WA_NdadaRs, WAUSA_NfiltRs, verbose=TRUE)
#head(WA_Nmergers[[1]])

#African American
AA_CRCmergers = mergePairs(AA_CRCdadaFs, AAUSA_CRCfiltFs, AA_CRCdadaRs, AAUSA_CRCfiltRs, verbose=TRUE)
AA_Nmergers = mergePairs(AA_NdadaFs, AAUSA_NfiltFs, AA_NdadaRs, AAUSA_NfiltRs, verbose=TRUE)

#part 6: Construction of amplicon sequence variant (ASV) table
#White American
WA_CRCseqtab = makeSequenceTable(WA_CRCmergers)
WA_Nseqtab = makeSequenceTable(WA_Nmergers)

#African American
# construct an amplicon sequence variant (ASV) table
AA_CRCseqtab = makeSequenceTable(AA_CRCmergers)
AA_Nseqtab = makeSequenceTable(AA_Nmergers)

dim(WA_Nseqtab)
table(nchar(getSequences(WA_Nseqtab))) #Inspect distribution of sequence lengths

# Part 7: removing chimeras
#White American
WA_CRCseqtab.nochim = removeBimeraDenovo(WA_CRCseqtab, method="consensus", verbose=TRUE)
WA_Nseqtab.nochim = removeBimeraDenovo(WA_Nseqtab, method="consensus", verbose=TRUE)
table(nchar(getSequences(WA_Nseqtab.nochim)))
table(nchar(getSequences(WA_CRCseqtab.nochim)))
sum(WA_CRCseqtab.nochim)/sum(WA_CRCseqtab)
sum(WA_Nseqtab.nochim)/sum(WA_Nseqtab)

#African American
AA_CRCseqtab.nochim = removeBimeraDenovo(AA_CRCseqtab, method="consensus", verbose=TRUE)
AA_Nseqtab.nochim = removeBimeraDenovo(AA_Nseqtab, method="consensus", verbose=TRUE)
table(nchar(getSequences(AA_Nseqtab.nochim)))
table(nchar(getSequences(AA_CRCseqtab.nochim)))
sum(AA_CRCseqtab.nochim)/sum(AA_CRCseqtab)
sum(AA_Nseqtab.nochim)/sum(AA_Nseqtab)

#part 8: Tracking the reads
getN = function(x) sum(getUniques(x))
#White American
WA_CRCtrack = cbind(WAUSACRC_out, sapply(WA_CRCdadaFs, getN), sapply(WA_CRCdadaRs, getN), sapply(WA_CRCmergers, getN), rowSums(WA_CRCseqtab.nochim))
WA_Ntrack = cbind(WAUSAN_out, sapply(WA_NdadaFs, getN), sapply(WA_NdadaRs, getN), sapply(WA_Nmergers, getN), rowSums(WA_Nseqtab.nochim))

#African American
AA_CRCtrack = cbind(AAUSACRC_out, sapply(AA_CRCdadaFs, getN), sapply(AA_CRCdadaRs, getN), sapply(AA_CRCmergers, getN), rowSums(AA_CRCseqtab.nochim))
AA_Ntrack = cbind(AAUSAN_out, sapply(AA_NdadaFs, getN), sapply(AA_NdadaRs, getN), sapply(AA_Nmergers, getN), rowSums(AA_Nseqtab.nochim))

#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#White American
colnames(WA_CRCtrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(WA_CRCtrack) = WAUSA_sample.names
colnames(WA_Ntrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(WA_Ntrack) = WAUSAN_sample.names

#African American
colnames(AA_CRCtrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(AA_CRCtrack) = AAUSA_sample.names
colnames(AA_Ntrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(AA_Ntrack) = AAUSAN_sample.names

#part 9: assigning taxonomy
WA_CRCtaxa = assignTaxonomy(WA_CRCseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
WA_Ntaxa = assignTaxonomy(WA_Nseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#inspect the taxonomic assignments
WA_CRCtaxa.print = WA_CRCtaxa #Removing sequence row names for display only
rownames(WA_CRCtaxa.print) = NULL
head(WA_CRCtaxa.print)
colnames(WA_CRCtaxa)

#African American
AA_CRCtaxa = assignTaxonomy(AA_CRCseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
AA_Ntaxa = assignTaxonomy(AA_Nseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#inspect the taxonomic assignments
AA_CRCtaxa.print = AA_CRCtaxa
rownames(AA_CRCtaxa.print) = NULL
head(WA_CRCtaxa.print)
colnames(WA_CRCtaxa)

#save.image(file = "DADA2_USA.rData")
