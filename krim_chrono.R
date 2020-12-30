ls | grep XXX | xargs cp -t  ~/storage/XXX/
ls | grep kimeklis_K_ | xargs cp -t  ~/storage/krim_chrono/in

setwd("~/storage/Ant_cr/16S/")

require(dada2)
require(Biostrings)
require(DECIPHER)
require(seqinr)
require(data.table)
require(metagMisc)
require(tibble)

path = "in/"
path_trein_set = "~/storage/scripts_git_done/dada2_pipeline/data/silva_nr_v132_train_set.fa"
path_trein_set_species = "~/storage/scripts_git_done/dada2_pipeline/data/silva_species_assignment_v132.fa"
name_Brief = "krimcronoBrief.txt"
truncLen = "210,180"
maxEE = "2,5"
mult = TRUE
mlt = NULL

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
on.exit()
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,180), trimLeft=c(20,19), maxN=0, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, multithread=mult)
errF <- learnErrors(filtFs, multithread=mult)
errR <- learnErrors(filtRs, multithread=mult)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=mult, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=mult, pool=TRUE)
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

write.table(briefToSeq, file = name_Brief, sep = "\t")

dna.ref <- DNAStringSet(briefToSeq)
alignment <- AlignSeqs(DNAStringSet(dna.ref), anchor=NA,verbose=FALSE, processors = mlt)
writeXStringSet(alignment, file="align.fasta")

#assign with DECIPER
# # with 1000 bootstrap
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# ids <- IdTaxa(dna, trainingSet, bootstraps = 1000, strand="top", processors=NULL, verbose=FALSE) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxa.dec <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa.dec <- x$taxon[m]
#   taxa.dec[startsWith(taxa.dec, "unclassified_")] <- NA
#   taxa.dec
# }))
# colnames(taxa.dec) <- ranks; rownames(taxa.dec) <- getSequences(seqtab.nochim)

#with 100 bootstrap

ids.100 <- IdTaxa(dna, trainingSet, bootstraps = 100, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
taxa.dec.100 <- t(sapply(ids.100, function(x) {
  m <- match(ranks, x$rank)
  taxa.dec.100 <- x$taxon[m]
  taxa.dec.100[startsWith(taxa.dec.100, "unclassified_")] <- NA
  taxa.dec.100
}))

colnames(taxa.dec.100) <- ranks; rownames(taxa.dec.100) <- getSequences(seqtab.nochim)

#assign with assign taxonomy
taxa.dada2 <- assignTaxonomy(briefToSeq,path_trein_set , multithread=mult, minBoot = 80)


taxa.dada2.species <- assignSpecies(briefToSeq, path_trein_set_species) # maybe use not seqs but brief 


rownames(taxa.dada2.species) <- rownames(briefToSeq)
briefToSeq.df <- data.frame(briefToSeq)
rownames(taxa.dada2.species) <- rownames(briefToSeq.df)

rownames(taxa.dada2) <- rownames(taxa.dada2.species)
taxa.rdp <- cbind2(taxa.dada2, taxa.dada2.species[,2])
colnames(taxa.rdp)[7] <- "Species"
colnames(taxa.rdp) <- ranks

rownames(taxa.dec) <- rownames(taxa.dada2.species)
taxa.dec <- cbind2(taxa.dec[,1:6], taxa.dada2.species[,2])
colnames(taxa.dec)[7] <- "Species"

rownames(taxa.dec.100) <- rownames(taxa.dada2.species)
taxa.dec.100 <- cbind2(taxa.dec.100[,1:6], taxa.dada2.species[,2])
colnames(taxa.dec.100)[7] <- "Species"

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file = "track.tsv", sep= "\t", col.names = NA, quote=FALSE)

st.brief.t.df <- data.frame(t(st.brief))
write.table(st.brief.t.df, file = "otu_table.txt", sep= "\t", col.names = NA, quote=FALSE)

briefToSeq.ls <- as.list(briefToSeq.df[,c("briefToSeq")])
briefToSeq.names <- as.list(rownames(briefToSeq.df))
write.fasta( briefToSeq.ls, briefToSeq.names , "rep_seq.fasta", as.string = FALSE)

write.table(taxa.rdp, file = "taxa.rdp.txt", sep= "\t", col.names = NA, quote=FALSE)
write.table(taxa.dec, file = "taxa.dec.txt", sep= "\t", col.names = NA, quote=FALSE)
write.table(taxa.dec.100, file = "taxa.dec.100.txt", sep= "\t", col.names = NA, quote=FALSE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file = "track.tsv", sep= "\t", col.names = NA, quote=FALSE)

# find the similarities/dissimilarities in the present assigns

head(taxa.rdp)
head(taxa.dec.100)
head(taxa.dec)

class(taxa.rdp)

library(tidyverse)
taxa.rdp.t <- taxa.rdp %>% as_tibble(rownames = "id")
taxa.dec.100.t <- taxa.dec.100 %>% as_tibble(rownames = "id")
taxa.dec.t <- taxa.dec %>% as_tibble(rownames = "id")

# list with the assign info
# from tibble
taxa.rdp.t %>% map(~ sum(!is.na(.))/length(.))
# from shitty matrix
taxa.dec %>% as.tibble() %>% map(~ sum(!is.na(.))/length(.))

# find spechiphyc pattern in GlobalEnv, then create a list of dataframes
ls(pattern = "taxa\\..+\\.t$", envir = .GlobalEnv) %>% 
  mget(envir = .GlobalEnv) 


%>% 
  tibble()
  
  
  map(~ sum(!is.na(.))/length(.))

distinct(taxa.rdp.t$Kingdom, taxa.dec.100.t$domain)
taxa.dec.100 %>% as.tibble(rownames = "id")
taxa.dec %>% as.tibble(rownames = "id") 

plot(ids)
par(new=TRUE)
plot(ids.100, add=TRUE)

par(mfrow=c(1,2))
plot(ids)
plot(ids.100)

trainingSet$problemGroups

require(dada2)
require(Biostrings)
require(DECIPHER)
require(phyloseq)
require(seqinr)
require(data.table)
require(metagMisc)
require(tibble)

setwd("~/storage/krim_chrono/shit_data/")

filt.otu <-t(as.data.frame(fread("filtered_table.txt")))
first <-  filt.otu[1,]
filt.otu <- filt.otu[-c(1),]
colnames(filt.otu) <-  first
class(filt.otu) <- "numeric"
filt.otu.matrix <- as.matrix(filt.otu)
#head.col <- scan("head.txt", character(), quote = "")
#rownames(filt.otu.matrix) <- head.col
tree <- read_tree(treefile="tree.nwk")

mapp <- read.csv("krim_agro_full.csv" , header=TRUE, sep="\t")
map <- data.frame(row.names="ID", mapp)

taxa <- read.csv("taxa.rdp.txt" , header=TRUE, sep="\t", row.names = "X")
taxa <- as.matrix(taxa)
ps <- phyloseq(otu_table(filt.otu.matrix, taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(taxa),
               phy_tree(tree))
ps

saveRDS(ps, file = "~/storage/plates/ps.krim_chrono.rds")
