
#first run - with controls(its), pool - true, 2.5ee, 220-180



path = "in/"
path_trein_set = "data/silva_nr_v132_train_set.fa"
path_trein_set_species = "data/silva_species_assignment_v132.fa"
name_Brief = "nameBrief.txt"
truncLen = "220,180"
maxEE = "2,5"
mult = TRUE
mlt = NULL
  require(dada2)
  require(Biostrings)
  require(DECIPHER)
  require(phyloseq)
  require(seqinr)
  require(data.table)
  require(metagMisc)
  require(tibble)




  fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  on.exit()
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(185,185), trimLeft=20, maxN=0, maxEE=c(6,8), rm.phix=TRUE, compress=TRUE, multithread=mult)
  errF <- learnErrors(filtFs, multithread=mult)
  errR <- learnErrors(filtRs, multithread=mult)
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  dadaFs <- dada(derepFs, err=errF, multithread=mult, pool=FALSE)
  dadaRs <- dada(derepRs, err=errR, multithread=mult, pool=FALSE)
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
  
  write.table(taxa, file = "taxa.txt", sep= "\t", col.names = NA, quote=FALSE)
  
  #some SEPP bash tree magic
  mapp <- read.csv("agr_map.csv" , header=TRUE, sep="\t")
  map <- data.frame(row.names="ID", mapp)
  
  taxa <- read.csv("taxa.txt" , header=TRUE, sep="\t")
  taxa <- column_to_rownames(taxa, 'X')
  taxa <- as.matrix(taxa)

  
  filt.otu <- read.csv("filtered_table.txt" , header=TRUE, sep="\t")
  colnames(filt.otu)[1] <- "ID"
  filt.otu <- column_to_rownames(filt.otu, "ID")
  filt.otu <- filt.otu[c(naturalsort(colnames(filt.otu)))]
  colnames(filt.otu) <- rownames(map)
  # class(filt.otu) <- "numeric"
  
  filt.otu <- as.matrix(filt.otu)

  #head.col <- scan("head.txt", character(), quote = "")
  #rownames(filt.otu.matrix) <- head.col
  tree <- read_tree(treefile="tree.nwk")
  

  
  ps <- phyloseq(otu_table(t(filt.otu), taxa_are_rows=FALSE), 
                 sample_data(map), 
                 tax_table(taxa),
                 phy_tree(tree))
  
  source("~/storage/scripts_git_done/functions.R")
  ps.f <- ps
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
  
  # data postprocessing. Delete and combine samples
  ps.f <- prune_samples(sample_names(ps.f) != c("mal_8"), ps.f)
  ps.f <- prune_samples(sample_names(ps.f) != c("mal_62"), ps.f)
  ps.f <- prune_samples(sample_names(ps.f) != c("mal_33"), ps.f)
  ps.f  <- prune_taxa(taxa_sums(ps.f) > 0, ps.f)
  
  
  ps.f.nex <- prune_samples(sample_data(ps.f)$Type %in% c("dna","cdna"), ps.f)
  ps.f.nex  <- prune_taxa(taxa_sums(ps.f.nex) > 0, ps.f.nex)
  
  ps.f.e <- prune_samples(sample_data(ps.f)$Horizon %in% c("E"), ps.f)
  ps.f.e  <- prune_taxa(taxa_sums(ps.f.e) > 0, ps.f.e)
  
  ps.f.e.5 <- prune_samples(sample_data(ps.f.e)$Site %in% c("5"), ps.f.e)
  ps.f.e.5  <- prune_taxa(taxa_sums(ps.f.e.5) > 0, ps.f.e.5)
  
  ps.f.b <- prune_samples(sample_data(ps.f)$Horizon %in% c("B"), ps.f)
  ps.f.b  <- prune_taxa(taxa_sums(ps.f.b) > 0, ps.f.b)
  
  ps.f.b.nex <- prune_samples(sample_data(ps.f.b)$Type %in% c("dna","cdna"), ps.f.b)
  ps.f.b.nex  <- prune_taxa(taxa_sums(ps.f.b.nex) > 0, ps.f.b.nex)
  
  ps.f.e.3 <- prune_samples(sample_data(ps.f.e)$Site %in% c("3"), ps.f.e)
  ps.f.e.3  <- prune_taxa(taxa_sums(ps.f.e.3) > 0, ps.f.e.3)
  
  ps.f.nex <- prune_samples(sample_data(ps.f)$Type %in% c("dna","cdna"), ps.f)
  ps.f.nex  <- prune_taxa(taxa_sums(ps.f.nex) > 0, ps.f.nex)
  
  ps.f.b.d <- prune_samples(sample_data(ps.f.b)$Type != c("cdna"), ps.f.b)
  ps.f.b.d  <- prune_taxa(taxa_sums(ps.f.b.d) > 0, ps.f.b.d)
  
  ps.f.nex.rna <- prune_samples(sample_data(ps.f.nex)$Type %in% c("cdna"), ps.f.nex)
  ps.f.nex.rna  <- prune_taxa(taxa_sums(ps.f.nex.rna) > 0, ps.f.nex.rna)
  
  ps.f.nex.dna <- prune_samples(sample_data(ps.f.nex)$Type %in% c("dna"), ps.f.nex)
  ps.f.nex.dna  <- prune_taxa(taxa_sums(ps.f.nex.dna) > 0, ps.f.nex.dna)
  
  ps.f.nex <- prune_samples(sample_data(ps.f)$Type %in% c("dna","cdna"), ps.f)
  ps.f.nex  <- prune_taxa(taxa_sums(ps.f.nex) > 0, ps.f.nex)
  
  
  

  amp.ne <- phyloseq_to_amp(ps.f.ne)
  
  amp<- phyloseq_to_amp(ps.f)
  amp_heatmap(amp, group_by = "Horizon")
  amp
  amp.all <- phyloseq_to_amp_without_tree(ps.f)

  amp_heatmap(amp, group_by = "Repeats",tax_add = "Phylum", tax_show = 40, tax_aggregate = "Family") + theme_bw() + theme(text = element_text(size=15), legend.position = "none")
  amp_heatmap(amp, group_by = "Repeats", tax_show = 15, tax_aggregate = "Phylum") + theme_bw() + theme(text = element_text(size=17), legend.position = "none")

  

  ordination.b <- ordinate(ps.f, "PCoA", "bray")
  p = plot_ordination(ps.f, ordination.b, type="sample", color="Repeats", title="PCoA - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p <- p + stat_ellipse( type="norm", alpha=0.7)
p

# permanova. adonis2. run in the console. needd some nice wrapper
  permanova.16cell <- function(ps, dist = "bray"){
    require(phyloseq)
    require(vegan)
    dist <- phyloseq::distance(ps, dist)
    metadata <- as(sample_data(ps), "data.frame")
    ad <- adonis2(dist ~ Fraction, data = metadata, permutations = 10000)
    return(ad)
  }
  
  ps.f.nex@sam_data$Type
  
  
  
  
  ps.14_26 <- prune_samples(sample_data(ps.f)$Repeats %in% c("d.14", "d.26"), ps.f)
  ps.14_26 <- prune_taxa(taxa_sums(ps.14_26) > 0, ps.14_26) 
permanova.16cell(ps.14_26)
# Repeats   1  1.69784 0.84194 53.266 0.0023 **

delete_eu <- function(ps){
  badTaxa <- taxa_names(subset_taxa(ps, Kingdom=="Eukaryota"))
  ps <- pop_taxa(ps, badTaxa)
  return(ps)
}

ps.f <- delete_eu(ps)

library(phyloseq)
ps.merged <- merge_samples(ps.f, "Repeats", fun = "median")
ps.merged <- prune_taxa(taxa_sums(ps.merged) > 0, ps.merged) 
estimate_richness(ps.merged)

beta_custom_norm_NMDS <- function(ps, seed = 7888, normtype="vst", color="Type", shape="Fraction"){
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
  p = plot_ordination(ps, ordination.b, type="sample", color, shape, label = "Replicate", title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color, shape, title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color, shape, title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1, common.legend = TRUE, legend = "right", font.label = list(size = 12, face = "bold", color ="black"))
  
  return(p.all)
}

beta_custom_norm_NMDS_elli <- function(ps, seed = 7888, normtype="vst", color="Repeats"){
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
p1 <- beta_custom_norm_NMDS_elli(ps.f.e.3, seed =5544)
p2 <- beta_custom_norm_NMDS_elli(ps.f.e.5, seed =5544)
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.99, label.x = 0.934,  labels = c("Точка 3","Точка 5"), font.label = list(size = 14, face = "bold", color ="black"))
beta_custom_norm_NMDS_elli(ps.f.e.5, seed =5544)

p1 <- beta_custom_norm_NMDS_elli(ps.f.e.3, seed =678)
p2 <- beta_custom_norm_NMDS_elli(ps.f.b.nex, seed =678)
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.992, label.x = 0.9,  labels = c("горизонт E","горизонт B"), font.label = list(size = 14, face = "bold", color ="black"))
beta_custom_norm_NMDS_elli(ps.f.e.5, seed =5544)

p1 <- beta_custom_norm_NMDS_elli(ps.f.e.3, seed =678)
p2 <- beta_custom_norm_NMDS_elli(ps.f.e.5, seed =678)
p3 <- beta_custom_norm_NMDS_elli(ps.f.b.nex, seed =678)
ggarrange(p1, p2, p3, ncol = 1 , nrow = 3, label.y = 0.993, label.x = 0.87,  labels = c("Точка 3 гор. E","Точка 5 гор. E","Точка 5 гор. B"), font.label = list(size = 14, face = "bold", color ="black"))

beta_custom_norm_NMDS_elli(ps.f.b, seed =1233)


#reads and ricnress
plot_rich_reads <- function(ps.f){
  rish <- estimate_richness(ps.f, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(ps.f))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Repeats"] <- ps.f@sam_data$Repeats
  library(ggrepel)
  p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus, y=reads)) + geom_text_repel(aes(x=otus, y=reads,label=Repeats)) + theme_bw() 
  library(ggforce)
  p1 <- p1  + geom_mark_ellipse(aes(x = otus, y=reads, group = Repeats, label = Repeats))
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
plot_rich_reads(ps.f.b)
plot_rich_reads_samlenames(ps.f.e)


#alpha plots
ps.f.i <- ps.f.e.3
alpha.custom <- function(ps.f.i, arrga = "PD"){
  require(picante)
  pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)
  
  pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
  pd <- cbind(pd,ps.f.i@sam_data, pd1 )
  
  bmi <- levels(pd$Подписи)
  bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
  p1 <- ggviolin(pd, x = "Подписи", y = arrga,
                 add = "boxplot", fill = "Подписи") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
    scale_x_discrete(limits = c( "ДНК", "ДНК<250","ДНК>250","кДНК", "кДНК<250", "кДНК>250")) + theme(axis.title.x = element_blank(), legend.title = element_blank())
  return(p1)
}

p.alpha.Cdt.oo <- alpha.custom(ps.f.e.3, arrga = "Observed")
p.alpha.Cdt.sh <- alpha.custom(ps.f.e.3, arrga = "Shannon")
p.alpha.Cdt.is <- alpha.custom(ps.f.e.3, arrga = "InvSimpson")
p.alpha.Cdt.pd <- alpha.custom(ps.f.e.3, arrga = "PD")

p.alpha.wt.oo <- alpha.custom(ps.f.e.5, arrga = "Observed")
p.alpha.wt.sh <- alpha.custom(ps.f.e.5, arrga = "Shannon")
p.alpha.wt.is <- alpha.custom(ps.f.e.5, arrga = "InvSimpson")
p.alpha.wt.pd <- alpha.custom(ps.f.e.5, arrga = "PD")

p.alpha.4.oo <- alpha.custom(ps.f.b.nex, arrga = "Observed")
p.alpha.4.sh <- alpha.custom(ps.f.b.nex, arrga = "Shannon")
p.alpha.4.is <- alpha.custom(ps.f.b.nex, arrga = "InvSimpson")
p.alpha.4.pd <- alpha.custom(ps.f.b.nex, arrga = "PD")

p1 <- ggarrange(p.alpha.Cdt.oo, p.alpha.Cdt.sh, p.alpha.Cdt.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p2 <- ggarrange(p.alpha.wt.oo, p.alpha.wt.sh,p.alpha.wt.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p3 <- ggarrange(p.alpha.4.oo, p.alpha.4.sh,p.alpha.4.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
ggarrange(p1, p2, p3, ncol = 1 , nrow = 3, label.y = 0.9855, label.x = 0.27,  labels = c("Точка 3 гор. E","Точка 5 гор. E","Точка 5 гор. B"), font.label = list(size = 14, face = "bold", color ="black"))  

p.alpha.Cdt.oo <- alpha.custom(ps.f.nex, arrga = "Observed")
p.alpha.Cdt.sh <- alpha.custom(ps.f.nex, arrga = "Shannon")
p.alpha.Cdt.is <- alpha.custom(ps.f.nex, arrga = "InvSimpson")
p.alpha.Cdt.pd <- alpha.custom(ps.f.nex, arrga = "PD")

ggarrange(p.alpha.Cdt.oo, p.alpha.Cdt.sh, p.alpha.Cdt.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)

alpha.custom <- function(ps.f.i, arrga = "PD"){
  require(picante)
  pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)
  
  pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
  pd <- cbind(pd,ps.f.i@sam_data, pd1 )
  
  bmi <- levels(pd$Хуёписи)
  bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
  p1 <- ggviolin(pd, x = "Хуёписи", y = arrga,
                 add = "boxplot", fill = "Хуёписи") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
    scale_x_discrete(limits = c( "E5", "E3","B3")) + theme(axis.title.x = element_blank(), legend.title = element_blank())
  return(p1)
}

p.alpha.Cdt.oo <- alpha.custom(ps.f.nex.rna, arrga = "Observed")
p.alpha.Cdt.sh <- alpha.custom(ps.f.nex.rna, arrga = "Shannon")
p.alpha.Cdt.is <- alpha.custom(ps.f.nex.rna, arrga = "InvSimpson")
p.alpha.Cdt.pd <- alpha.custom(ps.f.nex.rna, arrga = "PD")

p.alpha.wt.oo <- alpha.custom(ps.f.nex.dna, arrga = "Observed")
p.alpha.wt.sh <- alpha.custom(ps.f.nex.dna, arrga = "Shannon")
p.alpha.wt.is <- alpha.custom(ps.f.nex.dna, arrga = "InvSimpson")
p.alpha.wt.pd <- alpha.custom(ps.f.nex.dna, arrga = "PD")

p1 <- ggarrange(p.alpha.Cdt.oo, p.alpha.Cdt.sh, p.alpha.Cdt.is,p.alpha.Cdt.pd , ncol = 4 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p2 <- ggarrange(p.alpha.wt.oo, p.alpha.wt.sh,p.alpha.Cdt.is,p.alpha.wt.pd , ncol = 4 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.9855, label.x = 0.40,  labels = c("кДНК","ДНК"), font.label = list(size = 14, face = "bold", color ="black"))  


Des.Al <- function(ps){
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

des.b.rna <- Des.Al(ps.f.b)

View(des.b.rna)  
  
  
Des.Tax = function(ps, Taxa){
  require(microbiomeSeq)
  ps <- taxa_level(ps, Taxa)
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
  nice <- cbind(res.df, as.data.frame(aggdata)[rownames(res.df),])
  return(nice)
}  

des.b.rna.phylum <- Des.Tax(ps.f.b, "Phylum")
View(des.b.rna.phylum)


diagdds = phyloseq_to_deseq2(ps.f.e.3, ~ Type)                  
diagdds = estimateSizeFactors(diagdds, type="poscounts")
diagdds = estimateDispersions(diagdds, fitType = "local") 
diagdds = DESeq(diagdds)
samp <-sample_data(ps.f.e.3)
dds.counts <- diagdds@assays@.xData$data$counts
dds.counts.df <- as.data.frame(dds.counts)
aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$Type), median))
colnames(aggdata) <- aggdata[1,]
aggdata <- aggdata[-1,]
res = results(diagdds)
res.df <- as.data.frame(res)
nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),]) 



resscss <- function(x)(2*x^7-0.1)
dat <- resscss(ani.cor)
p1 <- ggcorrplot(dat, type = "lower", outline.col = "white", show.legend = FALSE, hc.order = TRUE)

resscss_back <- function(x)((0.651836*(10*x + 1)^(1/7))*100)
p1$plot_env$label <- resscss_back(p1$plot_env$label)
p1 + geom_text(aes(label = round(p1$plot_env$label, 1)))


amp_heatmap_grouped <- function(ps){
  require(ampvis2)
  require(phyloseq)
  require(ggpubr)
  
  # split the original dataset into two then convert the goddamnned phyloseq object into ampvis2 asshanded class
  ps.Cd <- prune_samples(sample_data(ps)$Хуёписи %in% c("E3"), ps)
  ps.Cd <- prune_taxa(taxa_sums(ps.Cd) > 0, ps.Cd)  
  amp.Cd <- phyloseq_to_amp(ps.Cd)
  
  ps.wCd <- prune_samples(sample_data(ps)$Хуёписи %in% c("E5"), ps)
  ps.wCd <- prune_taxa(taxa_sums(ps.wCd) > 0, ps.wCd)  
  amp.wCd <- phyloseq_to_amp(ps.wCd)
  
  ps.wwCd <- prune_samples(sample_data(ps)$Хуёписи %in% c("B3"), ps)
  ps.wwCd <- prune_taxa(taxa_sums(ps.wwCd) > 0, ps.wwCd)  
  amp.wwCd <- phyloseq_to_amp(ps.wwCd)
  
  # plot some ampvised heatmaps
  p.heat.Cd <- amp_boxplot(amp.wCd,
                           group_by = "Подписи",
                           tax_show = 8,
                           tax_aggregate = "Phylum",
                           tax_class = "Proteobacteria") + labs(title = "E3") + theme_bw() + theme(text = element_text(size=14))
  
  p.heat.wCd <- amp_boxplot(amp.Cd,
                            group_by = "Подписи",
                            tax_show = 8,
                            tax_aggregate = "Phylum",
                            tax_class = "Proteobacteria") + labs(title = "E5") + theme_bw()  + theme(text = element_text(size=14))
  
  p.heat.wwCd <- amp_boxplot(amp.wwCd,
                            group_by = "Подписи",
                            tax_show = 8,
                            tax_aggregate = "Phylum",
                            tax_class = "Proteobacteria") + labs(title = "B3") + theme_bw()  + theme(text = element_text(size=14))
  
  
  p <- ggarrange(p.heat.Cd, p.heat.wCd,p.heat.wwCd, ncol = 3,label.x = 0.105, nrow = 1, common.legend = TRUE) + theme(legend.title = element_blank())
  return(p)
}

p1 <- amp_heatmap_grouped(ps.f.nex)
p2 <- amp_heatmap_grouped(ps.f.i.Cdt)
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.985, label.x = 0.35,  labels = c("SGE.wt","SGE.Cdt"), font.label = list(size = 14, face = "bold", color ="black"))                   

