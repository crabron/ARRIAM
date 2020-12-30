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

des.otus <- Des.Al(ps.f.5)
View(des.otus)


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
