# dada2 pipe
# function dont't work

dada2_pipe <- function( path = "in/", path_trein_set = "data/RefSeq-RDPv2_16S_species.fa", path_trein_set_species = "data/RefSeq-RDP_dada2_assignment_species.fa", name_Brief = "rdp.plates.brief.txt", truncLen = "220,180", maxEE = "2,5", mult = TRUE, mlt = NULL){
require(dada2)
require(Biostrings)
require(DECIPHER)
require(phyloseq)
require(seqinr)
require(data.table)
require(metagMisc)
require(tibble)

#       fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#       fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
#       sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#       on.exit()
#       filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#       filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
      
#       out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180), trimLeft=20, maxN=0, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, multithread=mult)
#       errF <- learnErrors(filtFs, multithread=mult)
#       errR <- learnErrors(filtRs, multithread=mult)
#       derepFs <- derepFastq(filtFs, verbose=TRUE)
#       derepRs <- derepFastq(filtRs, verbose=TRUE)
#       # Name the derep-class objects by the sample names
#       names(derepFs) <- sample.names
#       names(derepRs) <- sample.names
#       dadaFs <- dada(derepFs, err=errF, multithread=mult, pool=TRUE)
#       dadaRs <- dada(derepRs, err=errR, multithread=mult, pool=TRUE)
#       mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#       seqtab <- makeSequenceTable(mergers)
#       dim(seqtab)
#       table(nchar(getSequences(seqtab)))
#       getN <- function(x) sum(getUniques(x))
#       seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=TRUE)
      
      
briefToSeq <- colnames(seqtab.nochim)
names(briefToSeq) <- paste0("Seq", seq(ncol(seqtab.nochim))) 
st.brief <- seqtab.nochim
colnames(st.brief) <- names(briefToSeq) 

write.table(briefToSeq, file = name_Brief, sep = "\t")

#dna <- DNAStringSet(briefToSeq)
#alignment <- AlignSeqs(DNAStringSet(dna), anchor=NA,verbose=FALSE, processors = mlt)
#writeXStringSet(alignment, file="align.fasta")

taxa.dada2 <- assignTaxonomy(briefToSeq,path_trein_set , multithread=mult)
taxa.dada2.species <- assignSpecies(briefToSeq, path_trein_set_species) # maybe use not seqs but brief 
rownames(taxa.dada2.species) <- rownames(briefToSeq)
briefToSeq.df <- data.frame(briefToSeq)
rownames(taxa.dada2.species) <- rownames(briefToSeq.df)
rownames(taxa.dada2) <- rownames(taxa.dada2.species)
taxa <- cbind2(taxa.dada2, taxa.dada2.species[,2])
colnames(taxa)[7] <- "Species"

st.brief.t.df <- data.frame(t(st.brief))
write.table(st.brief.t.df, file = "otu_table.txt", sep= "\t", col.names = NA, quote=FALSE)

briefToSeq.ls <- as.list(briefToSeq.df[,c("briefToSeq")])
briefToSeq.names <- as.list(rownames(briefToSeq.df))
write.fasta( briefToSeq.ls, briefToSeq.names , "rep_seq.fasta", as.string = FALSE)

# not working
system("conda activate qiime2-2019.1
sed -i '1s/^/#OTU_ID/' otu_table.txt 
biom convert -i otu_table.txt -o otu_table.biom --to-hdf5
qiime tools import \
--input-path rep_seq.fasta \
--output-path rep_seq.qza \
--type 'FeatureData[Sequence]'
qiime fragment-insertion sepp --i-representative-sequences rep_seq.qza  --o-tree insertion-tree.qza  --o-placements data/insertion-placements.qza --p-threads 30
qiime tools import \
--input-path otu_table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path feature-table.qza
qiime fragment-insertion filter-features \
--i-table feature-table.qza \
--i-tree insertion-tree.qza \
--o-filtered-table filtered_table.qza \
--o-removed-table removed_table.qza
unzip -p filtered_table.qza */data/* > filtered_table.biom
biom convert -i  filtered_table.biom -o filtered_table.txt --to-tsv
sed -i '1d' filtered_table.txt
unzip -p insertion-tree.qza */data/* > tree.nwk
       ")

filt.otu <-t(as.data.frame(fread("filtered_table.txt")))
first <-  filt.otu[1,]
filt.otu <- filt.otu[-c(1),]
colnames(filt.otu) <-  first
class(filt.otu) <- "numeric"
filt.otu.matrix <- as.matrix(filt.otu)
#head.col <- scan("head.txt", character(), quote = "")
#rownames(filt.otu.matrix) <- head.col
tree <- read_tree(treefile="tree.nwk")

#       mapp <- read.csv(map_path , header=TRUE, sep="\t")
#       map <- data.frame(row.names="ID", mapp)

ps.rdp <- phyloseq(otu_table(filt.otu.matrix, taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(taxa),
               phy_tree(tree))
return(ps.rdp)
}


filt.otu <-t(as.data.frame(fread("filtered_table.txt")))
first <-  filt.otu[1,]
filt.otu <- filt.otu[-c(1),]
colnames(filt.otu) <-  first 
class(filt.otu) <- "numeric"
filt.otu.matrix <- as.matrix(filt.otu)
#head.col <- scan("head.txt", character(), quote = "")
#rownames(filt.otu.matrix) <- head.col
tree <- read_tree(treefile="tree.nwk")

 mapp <- read.csv("new_plates.map.csv" , header=TRUE, sep=",")
 map <- data.frame(row.names="ID", mapp)

ps.rdp <- phyloseq(otu_table(filt.otu.matrix, taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(taxa),
               phy_tree(tree))
               
brief.raw <- read.csv("rdp.plates.brief.txt", sep="\t")
rdp.brief <- as.character(brief.raw$x)
names(rdp.brief) <- rownames(brief.raw)

taxa.dada2 <- assignTaxonomy(rdp.brief,"data/RefSeq-RDPv2_16S_species.fa" , multithread=TRUE)
taxa.dada2.species <- assignSpecies(rdp.brief, "data/RefSeq-RDP_dada2_assignment_species.fa", n=10000)  
rownames(taxa.dada2.species) <- rownames(rdp.brief)
briefToSeq.df <- data.frame(rdp.brief)
rownames(taxa.dada2.species) <- rownames(briefToSeq.df)
rownames(taxa.dada2) <- rownames(taxa.dada2.species)
rdp.taxa <- cbind2(taxa.dada2, taxa.dada2.species[,2])
colnames(taxa)[7] <- "Species"

filt.otu <-t(as.data.frame(fread("filtered_table.txt")))
first <-  filt.otu[1,]
filt.otu <- filt.otu[-c(1),]
colnames(filt.otu) <-  first
class(filt.otu) <- "numeric"
filt.otu.matrix <- as.matrix(filt.otu)
#head.col <- scan("head.txt", character(), quote = "")
#rownames(filt.otu.matrix) <- head.col
tree <- read_tree(treefile="tree.nwk")

mapp <- read.csv("new_plates.map.csv" , header=TRUE, sep=",")
map <- data.frame(row.names="ID", mapp)

ps <- phyloseq(otu_table(filt.otu.matrix, taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(taxa),
               phy_tree(tree))

# tree by IQtree
# not presented here


#align. DECIPHER. 
alignment <- AlignSeqs(DNAStringSet(Al.dna), anchor=NA,verbose=FALSE, processors = 4)
writeXStringSet(alignment, file="align.fasta")

#tree by iqtree
iqtree -s align.rdp.fasta -nt AUTO -ntmax 40 -bb 1000 -bnni 


      
# processing data
# popping fucking silt fraction. Soil science != science. Separate soils
ps.s.nosilt <- prune_samples(sample_data(ps.s.f)$Type != c("silt"), ps.s.f)
ps.s.black <- prune_samples(sample_data(ps.s.nosilt)$Soil %in% c("black soil"), ps.s.nosilt)
ps.s.pod <- prune_samples(sample_data(ps.s.nosilt)$Soil %in% c("sod-podzolic"), ps.s.nosilt)

#clear empty taxas
ps.s.nosilt <- prune_taxa(taxa_sums(ps.s.nosilt) > 0, ps.s.nosilt)
ps.s.black <- prune_taxa(taxa_sums(ps.s.black) > 0, ps.s.black)
ps.s.pod <- prune_taxa(taxa_sums(ps.s.pod) > 0, ps.s.pod)

#heatmap construction in ampvis2
amp.s.nosilt <- phyloseq_to_amp(ps.s.nosilt)
amp.s.black <- phyloseq_to_amp(ps.s.black)
amp.s.pod <- phyloseq_to_amp(ps.s.pod)
p.venn.allnosilt <- amp_venn(amp.s.nosilt, group_by = "Type", normalise = TRUE, cut_a = 0.000001, cut_f=0.000001)
p.venn.black <- amp_venn(amp.s.black, group_by = "Type", normalise = TRUE, cut_a = 0.000001, cut_f=0.000001)
p.venn.pod <- amp_venn(amp.s.pod, group_by = "Type", normalise = TRUE, cut_a = 0.000001, cut_f=0.000001)
p.venn.all <- ggarrange(p.venn.allnosilt, p.venn.black, p.venn.pod, labels = c("all","black_soil", "podzol"), ncol = 3 , nrow = 1)
p.venn.all

# beta PERMANOVA(aka adonis2 vegan)


# how the data is organised? 

#how about the back-up?      

save.image("~/storage/r-base_files/new_plates.RData")  

# beta, the little bit nicer one
      
beta_custom_norm_NMDS_elli <- function(ps, seed = 7888, normtype="vst", color="Repeats"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  require(ggforce)
  
  # beta_NMDS <- function(){
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  if (normtype =="vst")
    pst <- varianceStabilizingTransformation(diagdds)
  if (normtype =="log") 
    pst <- rlogTransformation(diagdds)
  
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  
  ps.rand <- rarefy_even_depth(ps, rngseed = seed)
  
  #beta and ordination
  
  ordination.b <- ordinate(ps.varstab, "NMDS", "bray")

  
  #plotting
  p1 = plot_ordination(ps, ordination.b, type="sample", color="пул_нуклеиновых_кислот", shape="почва", title="NMDS - Bray", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 20)) + geom_point(size = 5) +
    geom_mark_ellipse(aes(group = Подписи))
  
  #merge by ggpubr
  
  return(p1)
}
ps.s.nosilt.new <- ps.s.nosilt
ps.s.nosilt.new@sam_data <- as.data.frame(read.csv("map.nosilt.txt", sep = "\t"), row.names="ID")
ps.s.nosilt.new@sam_data <- read.csv("map.nosilt.txt", sep = "\t", row.names="ID")
mapp.nosilt <- read.csv("map.nosilt.txt", sep = "\t", row.names="ID")
sample_data(ps.s.nosilt.new) <-  mapp.nosilt
ps.s.nosilt.new
beta_custom_norm_NMDS_elli(ps.s.nosilt.new)

# try get phylotypes that cDNA > DNA > extcDNA, then plot the on some tree OMG

#for des soils
ps.s.nosilt.black <- prune_samples(sample_data(ps.s.nosilt.new)$Soil %in% c("black soil"), ps.s.nosilt.new)
ps.s.nosilt.black  <- prune_taxa(taxa_sums(ps.s.nosilt.black ) > 0, ps.s.nosilt.black) 

ps.s.nosilt.black.cdnadna <- prune_samples(sample_data(ps.s.nosilt.black)$Type != c("extracellular"), ps.s.nosilt.black)
ps.s.nosilt.black.cdnadna  <- prune_taxa(taxa_sums(ps.s.nosilt.black.cdnadna ) > 0, ps.s.nosilt.black.cdnadna) 

ps.s.nosilt.black.dnaextc <- prune_samples(sample_data(ps.s.nosilt.black)$Type != c("cdna"), ps.s.nosilt.black)
ps.s.nosilt.black.dnaextc  <- prune_taxa(taxa_sums(ps.s.nosilt.black.dnaextc ) > 0, ps.s.nosilt.black.dnaextc) 

ps.s.nosilt.pod <- prune_samples(sample_data(ps.s.nosilt.new)$Soil %in% c("sod-podzolic"), ps.s.nosilt.new)
ps.s.nosilt.pod  <- prune_taxa(taxa_sums(ps.s.nosilt.pod  ) > 0, ps.s.nosilt.pod ) 

ps.s.nosilt.pod.cdnadna  <- prune_samples(sample_data(ps.s.nosilt.pod)$Type != c("extracellular"), ps.s.nosilt.pod)
ps.s.nosilt.pod.cdnadna  <- prune_taxa(taxa_sums(ps.s.nosilt.pod.cdnadna ) > 0, ps.s.nosilt.pod.cdnadna) 

ps.s.nosilt.pod.dnaextc  <- prune_samples(sample_data(ps.s.nosilt.pod)$Type != c("cdna"), ps.s.nosilt.pod)
ps.s.nosilt.pod.dnaextc  <- prune_taxa(taxa_sums(ps.s.nosilt.pod.dnaextc ) > 0, ps.s.nosilt.pod.dnaextc) 

ps.s.nosilt.noext<- prune_samples(sample_data(ps.s.nosilt.new)$Type != c("extracellular"), ps.s.nosilt.new)
ps.s.nosilt.noext  <- prune_taxa(taxa_sums(ps.s.nosilt.noext ) > 0, ps.s.nosilt.noext) 

Des.soil.w.simper <- function(ps){
  require(DESeq2)
  require(vegan)
  require(tibble)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Type)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$Type), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
}      

# cDNA/DNA cDNA prevalense in the black soil

des.black.cdnadna <-  Des.soil.w.simper(ps.s.nosilt.black.cdnadna)
View(des.black.cdnadna)

filter.des.out.black.cdnadna <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <- des[(des$log2FoldChange < 0),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}

#des.otus.n1.pr <- des.otus.n1[rownames(des.otus.chr.n1[(is.na(des.otus.chr.n1$padj) != TRUE),]),]
#des.otus.n1.pr <- des.otus.n1.pr[rownames(des.otus.pod.n1[(is.na(des.otus.pod.n1$padj) != TRUE),]),]

des.black.filtered <- filter.des.out.black.cdnadna(des.black.cdnadna, alpha = 0.1, baseMean = 15)
length(rownames(des.black.filtered))
View(des.black.filtered)

# cDNA/DNA cDNA prevalense in the podzol soil

des.pod.cdnadna <-  Des.soil.w.simper(ps.s.nosilt.pod.cdnadna)
View(des.black.cdnadna)

filter.des.out.black.cdnadna <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <- des[(des$log2FoldChange < 0),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}

des.pod.filtered <- filter.des.out.black.cdnadna(des.pod.cdnadna, alpha = 0.1, baseMean = 15)
length(rownames(des.pod.filtered))
View(des.pod.filtered)

# making the fancy listy of rownames 

# list with some changed phylotypes

list.cdnadna <- c(rownames(des.pod.filtered), rownames(des.black.filtered)) 
list.soils <- c(list.cdnadna, rownames(des.black.ex.filtered), rownames(des.pod.ex.filtered))



pop.taxa <- function(physeq, badTaxa){
  require(phyloseq)
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

library(tidyverse)
library(phyloseq)
ps.cdnadna <- prune_taxa(list.cdnadna, ps.s.nosilt.noext )
ps.cdnadna.ex <-  prune_taxa(list.soils, ps.s.nosilt.noext )
ps.cdnadna
ps.cdnadna.ex

# list for the column taxa with the soils types occurence
list.soil <- ps.cdnadna@tax_table@.Data %>% as_tibble(rownames = "id") %>% 
  mutate(тип_почвы = case_when(id %in% intersect(rownames(des.black.filtered),rownames(des.pod.filtered)) ~ 'обе',
  id %in% rownames(des.pod.filtered) ~ 'подзол',
  id %in% rownames(des.black.filtered) ~ 'чернозём')) %>% 
  pull(тип_почвы)


taxa.cdnadna["тип_почвы"] <- list.soil 
taxa.cdnadna <-  ps.cdnadna@tax_table@.Data
taxa.cdnadna <- cbind(taxa.cdnadna, list.soil)
colnames(taxa.cdnadna)[8] <- "тип_почвы"
tax_table(ps.cdnadna) <- taxa.cdnadna

write.table(ps.cdnadna@tax_table, file = "taxa.cdnadna.txt", sep= "\t", col.names = NA, quote=FALSE)
ps.s.nosilt.new@sam_data <- read.csv("map.nosilt.txt", sep = "\t", row.names="ID")
mapp.nosilt <- read.csv("map.nosilt.txt", sep = "\t", row.names="ID")

ps <- ps.cdnadna
ps
require(phyloseq)
require(ggplot2)
require(ggpubr)
require(DESeq2)
require(ggforce)
library(ggrepel)
  
  # beta_NMDS <- function(){
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  if (normtype =="vst")
    pst <- varianceStabilizingTransformation(diagdds)
  if (normtype =="log") 
    pst <- rlogTransformation(diagdds)
  
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  

  pslog <- transform_sample_counts(ps, function(x) log(1 + x)) 
  ordination.b <- ordinate(pslog, "MDS", "bray")
  pslog@otu_table
  pslog@tax_table
  
  #plotting
  plot_ordination(ps, ordination.b, type="taxa",  color="Phylum", title="MDS - Bray", shape = "тип_почвы",
    axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 20)) + geom_point(size = 5) + 
    geom_text_repel(aes(label = Genus), size = 5, nudge_y = - 0.015)

beta_custom_norm_NMDS_species(ps.cdnadna)


# DeSeq2 for black soils
des.black.dnaextc <-  Des.soil.w.simper(ps.s.nosilt.black.dnaextc)

filter.des.out.black.dnaextc <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <- des[(des$log2FoldChange > 0),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}


des.black.ex.filtered <- filter.des.out.black.dnaextc(des.black.dnaextc, alpha = 0.1, baseMean = 50)
length(rownames(des.black.ex.filtered))
View(des.black.ex.filtered)

# for pod soil

des.black.dnaextc <-  Des.soil.w.simper(ps.s.nosilt.pod.dnaextc)
des.pod.ex.filtered <- filter.des.out.black.dnaextc(des.black.dnaextc, alpha = 0.1, baseMean = 50)
length(rownames(des.black.ex.filtered))


save.image("~/storage/r-base_files/new_plates.RData")

# what adout tree? only for cDNA increasing ASVs

ps.f.pr.soils <- prune_taxa(rownames(des.otus.n1.soils.pr), ps.f) 

ps.f.merged.soil <- merge_samples(ps.f.pr.soils, "Soil")
melted.soils <- psmelt(ps.f.merged.soil)
adundance <- melted.soils$Abundance / sum(melted.soils$Abundance)*100
sigtab = res[(res$padj < alpha), ]

phyloseq <-  ps.cdnadna

tree <- phyloseq@phy_tree
taxa.pruned <- as.data.frame(phyloseq@tax_table@.Data)
taxa.pruned <- taxa.pruned %>%  mutate_all(as.character)
taxa.pruned$number <- seq.int(from = nrow(taxa.pruned), to = 1)
taxa.pruned$taxa <- ifelse(is.na(taxa.pruned$Genus), taxa.pruned$Family, taxa.pruned$Genus)
taxa.pruned[taxa.pruned == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
taxa.pruned$taxa2 <- ifelse(is.na(taxa.pruned$Species), with(taxa.pruned, paste0(taxa)), with(taxa.pruned, paste0(taxa, " ", Species )))
taxa.pruned$taxa3 <- ifelse(taxa.pruned$Phylum == "Proteobacteria", with(taxa.pruned, paste0(taxa.pruned$number, ".", taxa2, " // ", Class)), with(taxa.pruned, paste0(taxa.pruned$number, ".", taxa2, " // ", Phylum)))
class(taxa.pruned$Kingdom)
tree$tip.label <- taxa.pruned$taxa3
p <- ggtree(tree, ladderize = F) + geom_tiplab(mapping = aes(), align=TRUE, linesize=.5) + xlim(NA, 4)
p

t.pod.cdna <-   as_tibble(des.pod.cdnadna, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(pod = log2FoldChange)

t.black.cdna <- as_tibble(des.black.cdnadna, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(black = log2FoldChange)




t.cdna <- full_join(t.pod.cdna, t.black.cdna) %>% mutate(pod = 0 - pod) %>% mutate(black = 0 - black)

all_equal(t.cdna$id, phyloseq@phy_tree$tip.label)

d.cdna <- as.data.frame(t.cdna) %>% column_to_rownames("id")
d.cdna.tips <- d.cdna
rownames(d.cdna.tips) <- tree$tip.label
d.cdna

p
p2 <- gheatmap(p, d.cdna.tips, offset=1, width=0.6, low = "white", high = "red", legend_title = "кДНК") 
p2

# adding exttracellular increasing ASVs to tree

phyloseq <-  ps.cdnadna.ex
taxa.pruned

tree <- phyloseq@phy_tree
taxa.pruned <- as.data.frame(phyloseq@tax_table@.Data)
taxa.pruned <- taxa.pruned %>%  mutate_all(as.character)
taxa.pruned$number <- seq.int(from = nrow(taxa.pruned), to = 1)
taxa.pruned$taxa <- ifelse(is.na(taxa.pruned$Genus), taxa.pruned$Family, taxa.pruned$Genus)
taxa.pruned[taxa.pruned == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
taxa.pruned$taxa2 <- ifelse(is.na(taxa.pruned$Species), with(taxa.pruned, paste0(taxa)), with(taxa.pruned, paste0(taxa, " ", Species )))
taxa.pruned$taxa3 <- ifelse(taxa.pruned$Phylum == "Proteobacteria", with(taxa.pruned, paste0(taxa.pruned$number, ".", taxa2, " // ", Class)), with(taxa.pruned, paste0(taxa.pruned$number, ".", taxa2, " // ", Phylum)))
taxa.pruned$taxa3 <- ifelse(is.na(taxa.pruned$taxa3),with(taxa.pruned, paste0(taxa.pruned$number, ".NA // NA")), taxa.pruned$taxa3)
class(taxa.pruned$Kingdom)
tree$tip.label <- taxa.pruned$taxa3
p <- ggtree(tree, ladderize = F) + geom_tiplab(mapping = aes(), align=TRUE, linesize=.5) + xlim(NA, 4)
p


# create tibbles for heatmaps
t.pod.cdna <-   as_tibble(des.pod.cdnadna, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(cdna_pod = log2FoldChange)

t.black.cdna <- as_tibble(des.black.cdnadna, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(cdna_black = log2FoldChange)
t.black.cdna
t.cdna <- full_join(t.pod.cdna, t.black.cdna) %>% mutate(cdna_pod = 0 - cdna_pod) %>% mutate(cdna_black = 0 - cdna_black)

t.black.ex <- as_tibble(des.black.ex.filtered, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(ex_black = log2FoldChange)

t.pod.ex <- as_tibble(des.pod.ex.filtered, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(ex_pod = log2FoldChange)

t.ex <- full_join(t.black.ex, t.pod.ex)
t.ex

t.cdna
t.ex
t.all <- full_join(t.cdna, t.ex)
t.all
d.all <- as.data.frame(t.all) %>% 
  column_to_rownames("id")
d.all.tips <- d.all  
tree$tip.label
rownames(d.all.tips) <- tree$tip.label

d.cdna.78 
  
d.cdna.78 <- d.all.tips %>% select(cdna_pod, cdna_black) %>% rename(дер.подзол = cdna_pod, чернозём = cdna_black)
d.ex.78 <- d.all.tips %>% select(ex_pod, ex_black) %>% rename(дер.подзол = ex_pod, чернозём = ex_black)

library(ggnewscale)
p2 <- p1 + new_scale_fill()
gheatmap(p2, df2, offset=15, width=.1,
         colnames_angle=90, colnames_offset_y = .25) +
  scale_fill_viridis_c(option="A", name="continuous\nvalue")



p2 <- gheatmap(p, d.cdna.78, offset=1.3, width=0.55, legend_title = "кДНК") 
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p2, d.ex.78, offset=2.3, width=0.55, legend_title = "экстрДНК") + scale_colour_gradientn(colours = c("blue", "white", "red"),na.value = "grey50", guide = "colourbar",aesthetics = "fill")
p3

# add base mean
# transform 
varstab <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats) # вставить своё                 
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)

  
  pst.dimmed <- t(as.matrix(assay(pst))) 
  # pst.dimmed[pst.dimmed < 0.0] <- 0.0 # если нужна ординация, то нужно т.к. нули
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  return(ps.varstab)
}

ps.s.nosilt.pod.varstab <- varstab(ps.s.nosilt.pod)
ps.s.nosilt.black.varstab <- varstab(ps.s.nosilt.black)
ps.s.nosilt.varstab <- varstab(ps.s.nosilt)



as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select_if(is.numeric) %>% 
  replace(is.na(.), 0) %>% 
  mutate(res = rowMeans(.)) %>% 
  pull(res) ->
  list.mean

as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select_if(is.numeric) %>% 
  replace(is.na(.), 0) %>% 
  rowwise() %>% 
  mutate(res=mean(c(.))) %>% 
  pull(res)
as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" )

#same for relative abundance
otus.pod.rel <- t(apply(otu_table(ps.s.nosilt.pod), 1, function(x) x / sum(x)))
otus.black.rel <- t(apply(otu_table(ps.s.nosilt.black), 1, function(x) x / sum(x)))

as_tibble(t(otus.pod.rel ), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select(id) -> id.pod
as_tibble(t(otus.pod.rel), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select_if(is.numeric) %>% 
  replace(is.na(.), 0) %>% 
  mutate(podMean=rowMeans(.)) %>% 
  select(podMean) -> res.pod
res.podMean.rel <- cbind2(id.pod,res.pod)
# for black
as_tibble(t(otus.black.rel), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select(id) -> id.black
as_tibble(t(otus.black.rel), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select_if(is.numeric) %>% 
  replace(is.na(.), 0) %>% 
  mutate(blackMean=rowMeans(.)) %>% 
  select(blackMean) -> res.black
res.blackMean.rel <- cbind2(id.black,res.black)

res.blackMean.rel
res.podMean.rel
all.resMean.rel <- full_join(res.podMean.rel, res.blackMean.rel)
all.resMean.rel %>% replace(is.na(.), 0) -> all.resMean.rel
all.resMean.rel$blackMean <- 0 - all.resMean.rel$blackMean

first.tree <- phyloseq@phy_tree
ids.tree <-  data.frame(id = first.tree$tip.label, id.taxa = tree$tip.label)
all.resMean.rel.melted
some.join <- full_join(ids.tree, all.resMean.rel.melted)
some.join <- subset(some.join, select = -c(id))
some.join
library(reshape2)
library(ggstance)
all.resMean.rel.melted <- melt(all.resMean.rel, id=c("id"))

some.join$value <- some.join$value*100
as.tibble(some.join)
p4 <- facet_plot(p3, panel = "дерново-подзолистая / чернозём", data = some.join,  geom = geom_barh, 
                 mapping = aes(value), color = "salmon", fill = "white",
                 stat='identity') 
all.resMean.rel.melted

library(gtable)
library(grid)
gt = ggplot_gtable(ggplot_build(p4))
gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.5*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
grid.draw(gt) # plot with grid draw


p3
p4


# selfmade varstab function - rowMean - not the optymal tidy solution
# for pod
as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select(id) -> id.pod
as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select_if(is.numeric) %>% 
  replace(is.na(.), 0) %>% 
  mutate(podMean=rowMeans(.)) %>% 
  select(podMean) -> res.pod
res.podMean <- cbind2(id.pod,res.pod)
# for black
as_tibble(t(ps.s.nosilt.black.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select(id) -> id.black
as_tibble(t(ps.s.nosilt.black.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  select_if(is.numeric) %>% 
  replace(is.na(.), 0) %>% 
  mutate(blackMean=rowMeans(.)) %>% 
  select(blackMean) -> res.black
res.blackMean <- cbind2(id.black,res.black)

res.blackMean
res.podMean
all.resMean <- full_join(res.podMean, res.blackMean)
all.resMean

p3


as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% pull(id) %>% length()



kimeklis.kdna.P1:zverev.mbm.plg.78

list.mean

as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>%
  filter(id %in% phyloseq@phy_tree$tip.label) %>%
  mutate(res = rowMeans(select_if(is.numeric)), na.rm = TRUE)


as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label)


 data.frame(list.sum, row.names = )
 
as_tibble(t(ps.s.nosilt.pod.varstab@otu_table@.Data), rownames = "id" ) %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) 


p4 <- facet_plot(p3, panel = '', data = d.cdna, geom = geom_barh, 
                 mapping = aes(pod, fill = black),
                 stat='identity', width = 0.8) 

# DivNet

library(DivNet)
divnet.res <-  divnet(ps.s.nosilt.pod, ncores = 30)
divnet.res

# not working at all

rownames(t.all.tips) <- tree$tip.label
tree$tip.label
d.cdna <- as.data.frame(t.cdna) %>% column_to_rownames("id")
d.cdna.tips <- d.cdna
rownames(d.cdna.tips) <- tree$tip.label

t.black.ex <- as_tibble(des.black.ex.filtered, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(black = log2FoldChange)

t.pod.ex <- as_tibble(des.pod.ex.filtered, rownames = "id") %>% 
  filter(id %in% phyloseq@phy_tree$tip.label) %>% 
  select(id, log2FoldChange) %>%  
  rename(pod = log2FoldChange)

t.ex <- full_join(t.black.ex, t.pod.ex)
d.ex <- as.data.frame(t.ex) %>% column_to_rownames("id")
d.ex.tips <- d.ex
rownames(d.ex.tips) <- tree$tip.label






p2 <- facet_plot(p, panel = '', data = d.cdna, geom = geom_barh, 
                 mapping = aes(pod, fill = black),
                 stat='identity', width = 0.8) 


pod = des.pod.cdnadna$log2FoldChange, chr = des.black.cdnadna$log2FoldChange

d3[d3<0] <- 0
pod.per <- d3$pod / (d3$pod + d3$chr)
chr.per <- d3$chr / (d3$pod + d3$chr) 
d3.norm <- data.frame(id=tree$tip.label, pod = pod.per, chr = chr.per)
df <- melt(d3.norm, id=c("id"))
df$value <- df$value*4
p2 <- facet_plot(p, panel = '', data = df, geom = geom_barh, 
                 mapping = aes(value, fill = variable),
                 stat='identity', width = 0.8) 

p2

library(grid)
library(gtable)

gt = ggplot_gtable(ggplot_build(p2))
gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.6*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
grid.draw(gt) # plot with grid draw

d3 <- data.frame(id=tree$tip.label, pod = des.otus.pod.n1.pr$log2FoldChange, chr = des.otus.chr.n1.pr$log2FoldChange)
d3$id <- rev(d3$id)
d3[d3<0] <- 0
d3$chr <- 0 - d3$chr
View(d3)
df <- melt(d3, id=c("id"))
df$value <- df$value/3
p2 <- facet_plot(p, panel = "дерново-подзолистая / чернозём", data = df, geom = geom_barh, 
                 mapping = aes(value, fill = variable),
                 stat='identity', width = 0.8) 

p2

library(grid)
library(gtable)

gt = ggplot_gtable(ggplot_build(p2))
gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.6*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
grid.draw(gt) # plot with grid draw

