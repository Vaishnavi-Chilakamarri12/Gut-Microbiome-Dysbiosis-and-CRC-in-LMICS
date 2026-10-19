# Install necessary packages
install.packages("BiocManager")
BiocManager::install("dada2", ask = FALSE)
install.packages("Rcpp")
library(Rcpp)
library(dada2)
library(ShortRead)
library(stats)

#Define paths for kenya data
path_CRC_ken = "D:/summer research MSU/SRA_data/Kenya/Ken_CRC"
path_N_ken = "D:/summer research MSU/SRA_data/Kenya/Ken_HC"
#list.files(path_N_ken)

#Sort and list FASTQ files
fnFs_CRC_ken = sort(list.files(path_CRC_ken, pattern="_1.fastq", full.names = TRUE))
fnFs_N_ken = sort(list.files(path_N_ken, pattern="_1.fastq", full.names = TRUE))
fnRs_CRC_ken = sort(list.files(path_CRC_ken, pattern="_2.fastq", full.names = TRUE))
fnRs_N_ken = sort(list.files(path_N_ken, pattern="_2.fastq", full.names = TRUE))

#Create a Vector to extract sample names up to the first underscore
KenCRC_sample.names = sapply(strsplit(basename(fnFs_CRC_ken), "_"), `[`, 1)
KenN_sample.names = sapply(strsplit(basename(fnFs_N_ken), "_"), `[`, 1)
KenCRC_sample.names = sapply(strsplit(basename(fnRs_CRC_ken), "_"), `[`, 1)
KenN_sample.names = sapply(strsplit(basename(fnRs_N_ken), "_"), `[`, 1)
fnFs_CRC_ken[1:7]

#Inspect read quality profiles
#plotQualityProfile(fnFs_CRC_ken[1:7])
#plotQualityProfile(fnRs_CRC_ken[1:7])
#plotQualityProfile(fnFs_N_ken[1:7])
#plotQualityProfile(fnRs_N_ken[1:7])

#part 2: Filter and trim data
ken_CRCfiltFs = file.path(path_CRC_ken, "filtered", paste0(KenCRC_sample.names, "_kenCRC_F_filt.fastq.gz"))
ken_CRCfiltRs = file.path(path_CRC_ken, "filtered", paste0(KenCRC_sample.names, "_kenCRC_R_filt.fastq.gz"))
ken_NfiltFs = file.path(path_N_ken, "filtered", paste0(KenN_sample.names, "_KenN_F_filt.fastq.gz"))
ken_NfiltRs = file.path(path_N_ken, "filtered", paste0(KenN_sample.names, "_KenN_R_filt.fastq.gz"))

names(ken_CRCfiltFs) = KenCRC_sample.names
names(ken_CRCfiltRs) = KenCRC_sample.names
names(ken_NfiltFs) = KenN_sample.names
names(ken_NfiltRs) = KenN_sample.names

kenCRC_out = filterAndTrim(fnFs_CRC_ken, ken_CRCfiltFs, fnRs_CRC_ken, ken_CRCfiltRs, 
                          truncLen = c(240, 150), maxN = 0, maxEE=c(2,5), truncQ=2, 
                          rm.phix = TRUE, compress = TRUE, multithread = FALSE)
KenN_out = filterAndTrim(fnFs_N_ken, ken_NfiltFs, fnRs_N_ken, ken_NfiltRs, 
                        truncLen = c(240, 150), maxN = 0, maxEE=c(2,5), truncQ=2, 
                        rm.phix = TRUE, compress = TRUE, multithread = FALSE)
#kenCRC_out
#KenN_out

#part 3: learn error rates
errF_kenCRC = learnErrors(ken_CRCfiltFs, multithread=TRUE)
errR_kenCRC = learnErrors(ken_CRCfiltRs, multithread=TRUE)
errF_KenN = learnErrors(ken_NfiltFs, multithread=TRUE)
errR_KenN = learnErrors(ken_NfiltRs, multithread=TRUE)

#Plot error rates for kenya samples
#plotErrors(errF_kenCRC, nominalQ=TRUE)
#plotErrors(errR_kenCRC, nominalQ = TRUE)
#plotErrors(errF_KenN, nominalQ=TRUE)
#plotErrors(errR_KenN, nominalQ = TRUE)
#black line-estimated error rates after convergence of the machine-learning algorithm
#red line- error rates expected under the nominal definition of the Q-score
#Here the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected.

#part 4: ASV Sample inference (used to filter the trimmed sequence data)
ken_CRCdadaFs = dada(ken_CRCfiltFs, err=errF_kenCRC, multithread=TRUE)
ken_CRCdadaRs = dada(ken_CRCfiltRs, err=errR_kenCRC, multithread=TRUE)
ken_NdadaFs = dada(ken_NfiltFs, err=errF_KenN, multithread=TRUE)
ken_NdadaRs = dada(ken_NfiltRs, err=errR_KenN, multithread=TRUE)

#part 5: Merge forward and reverse reads
ken_CRCmergers = mergePairs(ken_CRCdadaFs, ken_CRCfiltFs, ken_CRCdadaRs, ken_CRCfiltRs, verbose=TRUE)
ken_Nmergers = mergePairs(ken_NdadaFs, ken_NfiltFs, ken_NdadaRs, ken_NfiltRs, verbose=TRUE)

#part 6: Construction of amplicon sequence variant (ASV) table
ken_CRCseqtab = makeSequenceTable(ken_CRCmergers)
ken_Nseqtab = makeSequenceTable(ken_Nmergers)

# Part 7: removing chimeras
ken_CRCseqtab.nochim = removeBimeraDenovo(ken_CRCseqtab, method="consensus",verbose=TRUE)
table(nchar(getSequences(ken_CRCseqtab.nochim)))
sum(ken_CRCseqtab.nochim)/sum(ken_CRCseqtab)
ken_Nseqtab.nochim = removeBimeraDenovo(ken_Nseqtab, method="consensus", verbose=TRUE)
table(nchar(getSequences(ken_Nseqtab.nochim)))
sum(ken_Nseqtab.nochim)/sum(ken_Nseqtab)

#part 8: Tracking the reads
getN = function(x) sum(getUniques(x))
ken_CRCtrack = cbind(kenCRC_out, sapply(ken_CRCdadaFs, getN), sapply(ken_CRCdadaRs, getN), 
                  sapply(ken_CRCmergers, getN), rowSums(ken_CRCseqtab.nochim))
ken_Ntrack = cbind(KenN_out, sapply(ken_NdadaFs, getN), sapply(ken_NdadaRs, getN), 
                  sapply(ken_Nmergers, getN), rowSums(ken_Nseqtab.nochim))
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(ken_CRCtrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(ken_CRCtrack) = KenCRC_sample.names
colnames(ken_Ntrack) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(ken_Ntrack) = KenN_sample.names

#part 9: assigning taxonomy using SILVA database
ken_CRCtaxa = assignTaxonomy(ken_CRCseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
ken_Ntaxa = assignTaxonomy(ken_Nseqtab.nochim, "D:/summer research MSU/Tools/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#inspect the taxonomic assignments
ken_CRCtaxa.print = ken_CRCtaxa
rownames(ken_CRCtaxa.print) = NULL
head(ken_CRCtaxa.print)
colnames(ken_CRCtaxa)

#save.image(file = "DADA2_kenya.rData")
