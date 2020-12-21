library(dada2)
library(ShortRead)
library(Biostrings)
library(seqinr)

mult <- TRUE # TRUE
setwd("/home/gladkov/storage/cellulolit/its/")
path <- "raw"
# path <- "/home/gladkov/storage/cellulolit/its/"
path_trein_set <- "/home/gladkov/storage/pipez/data_1/sh_general_release_dynamic_01.12.2017.fasta"
FWD <- "GCATCGATGAAGAACGCAGC"  ## CHANGE ME to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"  ## CHANGE ME...

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))


allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, minLen = 100, multithread = TRUE, , maxEE = 2)

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/home/gladkov/.conda/envs/pandas/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
R3.flags <- paste("--minimum-lengt 10")
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags,R3.flags, "-n", 2,  # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files             
}

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 3), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

errF <- learnErrors(filtFs, multithread=mult)
errR <- learnErrors(filtRs, multithread=mult)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=mult)
dadaRs <- dada(derepRs, err=errR, multithread=mult)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
getN <- function(x) sum(getUniques(x))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=TRUE)

briefToSeq <- colnames(seqtab.nochim)
names(briefToSeq) <- paste0("Seq", seq(ncol(seqtab.nochim))) 
st.brief <- seqtab.nochim
colnames(st.brief) <- names(briefToSeq) 
st.brief.t.df <- data.frame(t(st.brief))
write.table(st.brief.t.df, file = "otu_table.txt", sep= "\t", col.names = NA, quote=FALSE)

write.table(briefToSeq, file = "briefToSeq.tsv", sep = "\t")

taxa <- assignTaxonomy(seqtab.nochim, path_trein_set, multithread = mult, tryRC = TRUE)
write.table(taxa, file = "taxa.txt", sep= "\t", col.names = NA, quote=FALSE)

briefToSeq.df <- data.frame(briefToSeq)
briefToSeq.ls <- as.list(briefToSeq.df[,c("briefToSeq")])
briefToSeq.names <- as.list(rownames(briefToSeq.df))
write.fasta( briefToSeq.ls, briefToSeq.names , "rep_seq.fasta", as.string = FALSE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
write.table(track, file = "track.tsv", sep= "\t", col.names = NA, quote=FALSE)