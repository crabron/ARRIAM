path_trein_set = "/home/gladkov/storage/scripts_git_done/dada2_pipeline/data/silva_nr_v132_train_set.fa"
path_trein_set_species = "/home/gladkov/storage/scripts_git_done/dada2_pipeline/data/silva_species_assignment_v132.fa"

#at first we have all base dada2 pipeline output, BUT we remake the original taxa that have been calc in DECIPHER package with original bd



# fasta to character vector
library("Biostrings")
library("tibble")

fastaFile <- readDNAStringSet("dna_rna_seqs.fna")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df <- column_to_rownames(df, "seq_name")
rep.char <- as.character(df$sequence)

#the taxa assign and colnames renaming 
taxa.dada2 <- assignTaxonomy( rep.char,path_trein_set , multithread=TRUE)
taxa.dada2.species <- assignSpecies( rep.char, path_trein_set_species)
taxa.dada2.df <- as.data.frame(taxa.dada2)
sum(is.na(taxa.dada2.df$Phylum))/length(rownames(taxa.dada2.df))*100

rownames(taxa.dada2.species) <- seq_name
briefToSeq.df <- data.frame(briefToSeq)
rownames(taxa.dada2.species) <- seq_name
rownames(taxa.dada2) <- rownames(taxa.dada2.species)
taxa <- cbind2(taxa.dada2, taxa.dada2.species[,2])
colnames(taxa)[7] <- "Species"

write.table(taxa, file = "taxa.txt", sep= "\t", col.names = NA, quote=FALSE)

#read the map, tree, otu_table from the outer files
library(phyloseq)
library(tibble)
mapp <- read.csv("metadata_dna_rna.txt", header=TRUE, sep="\t")
map <- data.frame(row.names="ID", mapp)

# taxa <- read.csv("taxa.txt" , header=TRUE, sep="\t")
# taxa <- column_to_rownames(taxa, 'X')
# taxa <- as.matrix(taxa)


filt.otu <- read.csv("dna_rna_count_table.txt" , header=TRUE, sep="\t")
filt.otu <- column_to_rownames(filt.otu, "X")
filt.otu <- as.matrix(filt.otu)

#head.col <- scan("head.txt", character(), quote = "")
#rownames(filt.otu.matrix) <- head.col
tree <- read_tree(treefile="tree.nwk")

ps <- phyloseq(otu_table(t(filt.otu), taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(taxa),
               phy_tree(tree))

ps
View(map)

pop_taxa <- function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

delete_mit_chl <- function(ps){
  badTaxa <- taxa_names(subset_taxa(ps, Order=="Chloroplast"))
  ps <- pop_taxa(ps, badTaxa)
  badTaxa <- taxa_names(subset_taxa(ps, Family=="Mitochondria"))
  ps <- pop_taxa(ps, badTaxa)
  return(ps)
}

ps.f <- delete_mit_chl(ps)

ps
ps.f

#reads and ricnress
plot_rich_reads <- function(ps.f){
  rish <- estimate_richness(ps.f, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(ps.f))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Repeats"] <- ps.f@sam_data$Repeats
  library(ggrepel)
  p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus, y=reads)) + geom_text_repel(aes(x=otus, y=reads,label=rownames(reads.summary))) + theme_bw()
  library(ggforce)
  p1 <- p1  + geom_mark_ellipse(aes(x = otus, y=reads, group = Repeats, label = Repeats), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  return(p1)
}

plot_rich_reads_samlenames <- function(ps.f){
  rish <- estimate_richness(ps.f, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(ps.f))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Repeats"] <- rownames(ps.f@sam_data)
  library(ggrepel)
  p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus, y=reads)) + geom_text_repel(aes(x=otus, y=reads,label=Repeats)) + theme_bw()
  return(p1)
}
plot_rich_reads(ps.f)
plot_rich_reads_samlenames(ps.f)


ps.f.dna <- prune_samples(sample_data(ps.f)$Pool %in% c("dna"), ps.f)
ps.f.dna <- prune_taxa(taxa_sums(ps.f.dna) > 0, ps.f.dna) 
ps.f.rna <- prune_samples(sample_data(ps.f)$Pool %in% c("rna"), ps.f)
ps.f.rna <- prune_taxa(taxa_sums(ps.f.rna) > 0,ps.f.rna)

ps.f.chr <- prune_samples(sample_data(ps.f)$Soil %in% c("chernozem"), ps.f)
ps.f.chr<- prune_taxa(taxa_sums(ps.f.chr) > 0, ps.f.chr) 
ps.f.pod <- prune_samples(sample_data(ps.f)$Soil %in% c("podzolic"), ps.f)
ps.f.pod <- prune_taxa(taxa_sums(ps.f.pod) > 0,ps.f.pod)


ps.f@sam_data$Soil

permanova.straw <- function(ps, dist = "bray"){
  require(phyloseq)
  require(vegan)
  dist <- phyloseq::distance(ps, dist)
  metadata <- as(sample_data(ps), "data.frame")
  ad <- adonis2(dist ~ Pool, data = metadata, permutations = 10000)
  return(ad)
}
permanova.straw(ps.f.chr, dist="bray")
permanova.straw(ps.f.pod, dist="bray")

permanova.straw(ps.f.chr, dist="unifrac")
permanova.straw(ps.f.pod, dist="unifrac")

permanova.straw(ps.f.chr, dist="wunifrac")
permanova.straw(ps.f.pod, dist="wunifrac")

permanova.straw <- function(ps, dist = "bray"){
  require(phyloseq)
  require(vegan)
  dist <- phyloseq::distance(ps, dist)
  metadata <- as(sample_data(ps), "data.frame")
  ad <- adonis2(dist ~ Duration, data = metadata, permutations = 10000)
  return(ad)
}

permanova.straw(ps.f.dna, dist="bray")
permanova.straw(ps.f.rna, dist="bray")

permanova.straw(ps.f.dna, dist="unifrac")
permanova.straw(ps.f.rna, dist="unifrac")

permanova.straw(ps.f.dna, dist="wunifrac")
permanova.straw(ps.f.rna, dist="wunifrac")

permanova.straw(ps.f, dist="bray")
permanova.straw(ps.f, dist="unifrac")
permanova.straw(ps.f, dist="wunifrac")

beta_custom_norm_NMDS_elli <- function(ps, seed = 4566, normtype="vst", color="Repeats"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  
  
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
  ordination.u <- ordinate(ps.rand, "NMDS", "unifrac")
  ordination.w <- ordinate(ps.varstab, "NMDS", "wunifrac")
  
  #plotting
  p1 = plot_ordination(ps, ordination.b, type="sample", color, title="NMDS - Bray", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Подписи, label = Подписи), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  p2 = plot_ordination(ps, ordination.u, type="sample", color, title="NMDS - unifrac", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Подписи, label = Подписи), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  p3 = plot_ordination(ps, ordination.w, type="sample", color, title="NMDS - wunifrac", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) + 
    geom_mark_ellipse(aes(group = Подписи, label = Подписи), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1, p2, p3, ncol = 3 , nrow = 1, common.legend = TRUE, legend = "none", font.label = list(size = 12, face = "bold", color ="black"))
  
  return(p.all)
}
library(ggpubr)
p1 <- beta_custom_norm_NMDS_elli(ps.f.pod)
p2 <- beta_custom_norm_NMDS_elli(ps.f.chr)
ggarrange(p1, p2, ncol = 1, nrow = 2, label.y = 0.99, label.x = 0.8,  labels = c("Дерново-подзолистая почва","Чернозем"), font.label = list(size = 14, face = "bold", color ="black"))

ps.f.5 <- prune_samples(sample_data(ps.f)$Duration %in% c("5"), ps.f)
ps.f.5 <- prune_taxa(taxa_sums(ps.f.5) > 0, ps.f.5) 

ps.f.chr.5 <- prune_samples(sample_data(ps.f.chr)$Duration %in% c("5"), ps.f.chr)
ps.f.chr.5 <- prune_taxa(taxa_sums(ps.f.5) > 0, ps.f.5) 

ps.f.pod.5 <- prune_samples(sample_data(ps.f.pod)$Duration %in% c("5"), ps.f.pod)
ps.f.pod.5 <- prune_taxa(taxa_sums(ps.f.pod.5) > 0, ps.f.pod.5) 

ps.f.n1 <- prune_samples(sample_data(ps.f)$Duration != c("1"), ps.f)
ps.f.n1 <- prune_taxa(taxa_sums(ps.f.n1) > 0, ps.f.n1) 

ps.f.chr.n1 <- prune_samples(sample_data(ps.f.chr)$Duration != c("1"), ps.f.chr)
ps.f.chr.n1 <- prune_taxa(taxa_sums(ps.f.chr.n1) > 0, ps.f.chr.n1) 

ps.f.pod.n1 <- prune_samples(sample_data(ps.f.pod)$Duration != c("1"), ps.f.pod)
ps.f.pod.n1 <- prune_taxa(taxa_sums(ps.f.pod.n1) > 0, ps.f.pod.n1) 



Des.Tax = function(ps, Taxa){
  require(DESeq2)
  require(phyloseq)
  require(microbiomeSeq)
  ps <- taxa_level(ps, Taxa)
  diagdds = phyloseq_to_deseq2(ps, ~ Pool)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Pool), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df, as.data.frame(aggdata)[rownames(res.df),])
  return(nice)
}  

des.phylo <- Des.Tax(ps.f.5, "Phylum")
View(des.phylo)

Des.Al <- function(ps){
  require(DESeq2)
  require(vegan)
  require(tibble)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Pool)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  otus.ps.vegan <- veganifyOTU(ps.varstab)
  metadata <- as(sample_data(ps.varstab), "data.frame")
  sim <- with(metadata, simper(otus.ps.vegan, Pool))
  simper.vst <- cbind(sim$dna_rna$species,sim$dna_rna$average)
  simper.vst <- as.data.frame(simper.vst)
  simper.vst <- column_to_rownames(simper.vst, "V1")
  
  otus.ps.vegan <- veganifyOTU(ps)
  metadata <- as(sample_data(ps), "data.frame")
  sim <- with(metadata, simper(otus.ps.vegan, Pool))
  simper <- cbind(sim$dna_rna$species,sim$dna_rna$average)
  simper <- as.data.frame(simper)
  simper <- column_to_rownames(simper, "V1")
  
  diagdds = phyloseq_to_deseq2(ps, ~ Pool)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$Pool), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),], simper[rownames(res.df),], simper.vst[rownames(res.df),])               
  return(nice)
}      

des.soils <-  
des.otus <- Des.Al(ps.f.5)
View(des.otus)

des.otus.pod.5 <- Des.Al(ps.f.pod.5)
des.otus.chr.5 <- Des.Al(ps.f.chr.5)

des.otus.n1 <- Des.Al(ps.f.n1)
des.otus.pod.n1 <- Des.Al(ps.f.pod.n1)
des.otus.chr.n1 <- Des.Al(ps.f.chr.n1)


otus.ps.vegan <- veganifyOTU(ps.f)
metadata <- as(sample_data(ps.f), "data.frame")
sim <- with(metadata, simper(otus.ps.vegan, Pool))
simper <- cbind(sim$s1903_In_s1903_In_Al$species,sim$s1903_In_s1903_In_Al$average)

otus.ps.vegan <- veganifyOTU(ps)
metadata <- as(sample_data(ps), "data.frame")
sim <- with(metadata, simper(otus.ps.vegan, Pool))
simper <- cbind(sim$dna_rna$species,sim$dna_rna$average)
simper <- as.data.frame(simper)

View(simper)

write.table(des.otus, file = "Des.5.txt", sep= "\t", col.names = NA, quote=FALSE)

alpha.custom <- function(ps.f.i, arrga = "PD"){
  require(picante)
  pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)
  
  pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
  pd <- cbind(pd,ps.f.i@sam_data, pd1 )
  
  bmi <- levels(as.factor(pd$Duration))
  bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
  p1 <- ggviolin(pd, x = "Duration", y = arrga,
                 add = "boxplot", fill = "Duration") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
    scale_x_discrete(limits = c("1", "5","10", "14")) + theme(axis.title.x = element_blank(), legend.title = element_blank())
  return(p1)
}

p.alpha.Cdt.oo <- alpha.custom(ps.f.dna, arrga = "Observed")
p.alpha.Cdt.sh <- alpha.custom(ps.f.dna, arrga = "Shannon")
p.alpha.Cdt.is <- alpha.custom(ps.f.dna, arrga = "InvSimpson")
p.alpha.Cdt.pd <- alpha.custom(ps.f.dna, arrga = "PD")

p.alpha.wt.oo <- alpha.custom(ps.f.rna, arrga = "Observed")
p.alpha.wt.sh <- alpha.custom(ps.f.rna, arrga = "Shannon")
p.alpha.wt.is <- alpha.custom(ps.f.rna, arrga = "InvSimpson")
p.alpha.wt.pd <- alpha.custom(ps.f.rna, arrga = "PD")

p1 <- ggarrange(p.alpha.Cdt.oo, p.alpha.Cdt.sh,p.alpha.Cdt.is, p.alpha.Cdt.pd , ncol = 4 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p2 <- ggarrange(p.alpha.wt.oo, p.alpha.wt.sh,p.alpha.wt.is, p.alpha.wt.pd , ncol = 4 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.9855, label.x = 0.32,  labels = c("ДНК","кДНК"), font.label = list(size = 14, face = "bold", color ="black"))  


require(picante)
pd <- pd(ps.f@otu_table@.Data, ps.f@phy_tree)

pd1 <- estimate_richness(ps.f, measures=c("Observed", "InvSimpson", "Shannon"))
pd <- cbind(pd,ps.f@sam_data, pd1 )

bmi <- levels(as.factor(pd$Duration))
bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])



library(tidyverse)
library(ggtree)

des.des.otus

#tree visualisation

ps.f.5.pr <- prune_taxa(taxa_sums(ps.f.5) > 1000, ps.f.5) 
ps.f.5.pr
tree.5.pr <- ps.f.5.pr@phy_tree
ggtree(tree.5.pr) + geom_tiplab()

plot_tree(tree.5.pr, color="SampleType", shape="Family", label.tips="Genus", size="Abundance")

require(ggplot2)
require(ggtree)

draw_phyloseq_ggtree <- function(phyloseq) {
  require(phyloseq)
  require(ggtree)
  library(scales)
  tree <- phyloseq@phy_tree
  p <- ggtree(tree, ladderize = F)
  # p <- p + geom_text(subset=.(!isTip), aes(label=label), hjust=-.2, size=4)
  
  dd <- psmelt(phyloseq)
  dd <- dd[dd$Abundance > 0, ]
  data <- merge(p$data, dd, by.x="label", by.y="OTU")
  
  spacing <- 0.02
  idx <- with(data, sapply(table(node)[unique(node)], function(i) 1:i)) %>% unlist
  hjust <- spacing * idx * max(data$x)
  data$xdodge <- data$x + hjust
  
  p <- p + geom_point(data=data, aes(x=xdodge, color=Pool,
                                     shape=Family, size=Abundance), na.rm=T) + 
    theme(legend.position="right") + scale_size_continuous(trans=log_trans(5))
  
  d2 <- unique(data[, c("x", "y", "Genus")])
  p + geom_text(data=d2, aes(label=Genus), hjust=-.3, na.rm=T, size=4)
}

draw_phyloseq_ggtree(ps.f.5.pr)


phyloseq <-  ps.f.pr
require(phyloseq)
require(ggtree)
library(scales)
library(ggplot2)
library(ggstance)
tree <- phyloseq@phy_tree
p <- ggtree(tree, ladderize = F)
# p <- p + geom_text(subset=.(!isTip), aes(label=label), hjust=-.2, size=4)
dd <- psmelt(phyloseq)
dd <- dd[dd$Abundance > 0, ]
data <- merge(p$data, dd, by.x="label", by.y="OTU")

spacing <- 0.02
data$Pool <- as.numeric(data$Pool)
hjust <- spacing * max(data$x) * data$Pool
data$xdodge <- data$x + hjust

d2 <- unique(data[, c("x", "y", "Genus")])
d4 <- unique(data[, c("x", "y", "Class")])
p <- p + geom_text(data=d2, aes(label=Genus), hjust=-.2, na.rm=T, size=4) 
p <- p + geom_text(data=d4, aes(label=Class), hjust=-.7, na.rm=T, size=4) 
p
p2 <- ggplot() + geom_point(data=data, aes(x=xdodge, color=Sample, size=Abundance), na.rm=T) + theme(legend.position="right") + scale_size_continuous(trans=log_trans(5))

des.otus.pruned <- des.otus[tree$tip.label,]
View(des.otus.pruned)

filter.des.out <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <-  des[(des$log2FoldChange > 0),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}

des.otus.n1.pr <- des.otus.n1[rownames(des.otus.chr.n1[(is.na(des.otus.chr.n1$padj) != TRUE),]),]
des.otus.n1.pr <- des.otus.n1.pr[rownames(des.otus.pod.n1[(is.na(des.otus.pod.n1$padj) != TRUE),]),]

des.otus.n1.pr <- filter.des.out(des.otus.n1.pr)
View(des.otus.chr.n1.pr)

des.otus.chr.n1.pr <- des.otus.chr.n1[rownames(des.otus.n1.pr),]
des.otus.pod.n1.pr <- des.otus.pod.n1[rownames(des.otus.n1.pr),]

ps.f.pr <- prune_taxa(rownames(des.otus.chr.n1.pr), ps.f) 


sigtab = res[(res$padj < alpha), ]


phyloseq <-  ps.f.pr
require(phyloseq)
require(ggtree)
library(scales)
library(ggplot2)
library(ggstance)
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


d3 <- data.frame(id=tree$tip.label, pod = des.otus.pod.n1.pr$log2FoldChange, chr = des.otus.chr.n1.pr$log2FoldChange)
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

library(ggpubr)
source("storage/scripts_git_done/functions.R")
amp <- phyloseq_to_amp(ps.f)
amp_heatmap(amp, group_by = "Pool",tax_add = "Phylum", tax_show = 40, tax_aggregate = "Family") + theme_bw() + theme(text = element_text(size=15), legend.position = "none")
amp_ph <- amp_heatmap(amp, group_by = "Pool", tax_show = 15, tax_aggregate = "Phylum") + theme_bw() + theme(text = element_text(size=17), legend.position = "none")
amp_class <- amp_heatmap(amp, group_by = "Pool", tax_show = 15, tax_aggregate = "Class") + theme_bw() + theme(text = element_text(size=17), legend.position = "none")
ggarrange(amp_ph, amp_class, ncol = 2 , nrow = 1)

dd <- psmelt(ps.f.pr)
sum(dd$Abundance)/1198446 * 100
dd <- psmelt(ps.f)
amp


#for des soils
Des.soil.w.simper <- function(ps){
  require(DESeq2)
  require(vegan)
  require(tibble)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Soil)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$Soil), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
}      

des.soils <-  Des.soil.w.simper(ps.f.n1)

filter.des.out <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}

des.otus.n1.pr <- des.otus.n1[rownames(des.otus.chr.n1[(is.na(des.otus.chr.n1$padj) != TRUE),]),]
des.otus.n1.pr <- des.otus.n1.pr[rownames(des.otus.pod.n1[(is.na(des.otus.pod.n1$padj) != TRUE),]),]

des.otus.n1.soils.pr <- filter.des.out(des.soils, alpha = 0.05, baseMean = 40)
length(rownames(des.otus.n1.soils.pr))

ps.f.pr.soils <- prune_taxa(rownames(des.otus.n1.soils.pr), ps.f) 

ps.f.merged.soil <- merge_samples(ps.f.pr.soils, "Soil")
melted.soils <- psmelt(ps.f.merged.soil)
adundance <- melted.soils$Abundance / sum(melted.soils$Abundance)*100
sigtab = res[(res$padj < alpha), ]

phyloseq <-  ps.f.pr.soils
require(phyloseq)
require(ggtree)
library(scales)
library(ggplot2)
library(ggstance)
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

df <- data.frame(id=tree$tip.label, value = des.otus.n1.soils.pr$log2FoldChange)
df$value <- df$value/7
p2 <- facet_plot(p, panel = "дерново-подзолистая / чернозём", data = df, geom = geom_barh, 
                 mapping = aes(value),
                 stat='identity', width = 0.8) 

# p2 <- facet_plot(p2, panel = "Abundance", data = df, geom = geom_barh, 
                 # mapping = aes(value),
                 # stat='identity', width = 0.8) 
p2
gt = ggplot_gtable(ggplot_build(p2))
gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.6*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
grid.draw(gt) # plot with grid draw


# some its

setwd("storage/straw_dna_rna/its/")
library(phyloseq)
library(tibble)
mapp <- read.csv("mapping_its_4.csv", header=TRUE, sep="\t")
map.its <- data.frame(row.names="ID", mapp)
taxa.its <- read.csv("taxa.txt" , header=TRUE, sep="\t")
brief.its <- read.csv("briefToSeq.tsv" , header=TRUE, sep="\t")
taxa.its <- as.data.frame(taxa.its)
brief.its <- as.data.frame(brief.its)
row.names(taxa.its) <-  row.names(brief.its)
taxa.its <- as.matrix(taxa.its)


otu.its <- read.csv("otu_table.txt" , header=TRUE, sep="\t")
otu.its <- column_to_rownames(otu.its, "X")
otu.its <- as.matrix(otu.its)

#head.col <- scan("head.txt", character(), quote = "")
#rownames(filt.otu.matrix) <- head.col

ps.its <- phyloseq(otu_table(t(otu.its), taxa_are_rows=FALSE), 
               sample_data(map.its), 
               tax_table(taxa.its))

ps.its <- prune_samples(sample_data(ps.its)$Work %in% c("yes"), ps.its)
ps.its <- prune_taxa(taxa_sums(ps.its) > 0, ps.its) 
amp.its <- phyloseq_to_amp_without_tree(ps.its)
amp.its
library(ampvis2)
order_x=naturalsort(levels(ps.its@sam_data$Подписи))
order_x_2 <- append(order_x, order_x)
amp_heatmap(amp.its, group_by = "Подписи",tax_add = "Phylum", tax_show = 40, tax_aggregate = "Genus", facet_by = "Почва") + theme_bw() + theme(text = element_text(size=15), legend.position = "none")
amp_ph <- amp_heatmap(amp, group_by = "Pool", tax_show = 15, tax_aggregate = "Phylum") + theme_bw() + theme(text = element_text(size=17), legend.position = "none")
amp_class <- amp_heatmap(amp, group_by = "Pool", tax_show = 15, tax_aggregate = "Class") + theme_bw() + theme(text = element_text(size=17), legend.position = "none")
amp.its
library(phyloseq)
ps.its@sam_data
beta_custom_norm_NMDS_elli_bray <- function(ps, normtype="vst", color="Подписи"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  require(ggforce)
  
  # beta_NMDS <- function(){
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Подписи)                  
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
  

  #beta and ordination
  
  ordination.b <- ordinate(ps.varstab, "PCoA", "bray")

  
  #plotting
  p1 = plot_ordination(ps, ordination.b, type="sample", color, title="PCoA - Bray", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Подписи, label = Подписи), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  p2 = plot_ordination(ps, ordination.b, type="sample", color, title="PCoA - Bray", 
                       axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Подписи, label = Подписи), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  

  #merge by ggpubr
  p.all <- ggarrange(p1, p2, ncol = 2 , nrow = 1, common.legend = TRUE, legend = "none", font.label = list(size = 12, face = "bold", color ="black"))
  return(p.all)
}
ps.its.chr <- prune_samples(sample_data(ps.its)$Soil %in% c("chernozem"), ps.its)
ps.its.chr<- prune_taxa(taxa_sums(ps.its.chr) > 0, ps.its.chr) 
ps.its.pod <- prune_samples(sample_data(ps.its)$Soil %in% c("soddy-podzolic"), ps.its)
ps.its.pod <- prune_taxa(taxa_sums(ps.its.pod) > 0,ps.its.pod)
library(ggpubr)
p1 <- beta_custom_norm_NMDS_elli_bray(ps.its.pod)
p2 <- beta_custom_norm_NMDS_elli_bray(ps.its.chr)
ggarrange(p1, p2, ncol = 1, nrow = 2, label.y = 0.99, label.x = 0.7,  labels = c("Дерново-подзолистая почва","Чернозем"), font.label = list(size = 14, face = "bold", color ="black"))

mapp_5 <- read.csv("~/storage/straw_dna_rna/its/mapping_its_5.csv", header=TRUE, sep="\t")
map.its.5 <- data.frame(row.names="ID", mapp_5)

amp.its
sam_data(ps.its) <- map.its.5 
source("~/storage/scripts_git_done/functions.R")
amp.n <- phyloseq_to_amp_without_tree(ps.its)
reshape2::melt(tree)

write.table(ps.s@otu_table, file = "otu.tsv", sep= "\t", col.names = NA, quote=FALSE)
write.table(ps.s@tax_table@.Data, file = "taxa.tsv", sep= "\t", col.names = NA, quote=FALSE)
write.table(ps.s@sam_data, file = "map.tsv", sep= "\t", col.names = NA, quote=FALSE)
ape::write.tree(ps.s@phy_tree, "tree.tree")