
#first run - with controls(its), pool - true, 2.5ee, 220-180



path = "in/"
path_trein_set = "data/silva_nr_v132_train_set.fa"
path_trein_set_species = "data/silva_species_assignment_v132.fa"
name_Brief = "cellulo.nameBrief.txt"
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
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,160), trimLeft=20, maxN=0, maxEE=c(2,3), rm.phix=TRUE, compress=TRUE, multithread=mult)
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
  briefToSeq <- read.csv("briefToSeq.tsv" , header=TRUE, sep="\t")
  taxa <- read.csv("taxa.txt" , header=TRUE, sep="\t")
  taxa <- column_to_rownames(taxa, "X")
  taxa <- as.matrix(taxa)
  rownames(taxa) <- rownames(briefToSeq)
  mapp <- read.csv("map.cel.16.csv" , header=TRUE, sep="\t")
  map <- data.frame(row.names="ID", mapp)
  
  
  filt.otu.f <- read.csv("otu_table.txt" , header=TRUE, sep="\t")
  colnames(filt.otu.f)[1] <- "ID"
  filt.otu <- column_to_rownames(filt.otu, "ID")
  filt.otu <- filt.otu[c(naturalsort(colnames(filt.otu)))]
  colnames(filt.otu) <- rownames(map)
  filt.otu <- filt.otu[-c(28, 27, 26, 25)]

  # class(filt.otu) <- "numeric"
  
  filt.otu <- as.matrix(filt.otu)

  #head.col <- scan("head.txt", character(), quote = "")
  #rownames(filt.otu.matrix) <- head.col

  

  
  ps <- phyloseq(otu_table(t(filt.otu), taxa_are_rows=FALSE), 
                 sample_data(map), 
                 tax_table(taxa))
  
  
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
  ps.f <- prune_samples(sample_data(ps.f)$Control %in% c("no"), ps.f)
  ps.f <- prune_taxa(taxa_sums(ps.f) > 0, ps.f)
  
  amp <- phyloseq_to_amp_without_tree(ps.f)
  amp.all <- phyloseq_to_amp_without_tree(ps.f)

  amp_heatmap(amp, group_by = "Repeats",tax_add = "Phylum", tax_show = 40, tax_aggregate = "Genus") + theme_bw() + theme(text = element_text(size=15), legend.position = "none")
  amp_heatmap(amp, group_by = "Repeats", tax_show = 15, tax_aggregate = "Phylum") + theme_bw() + theme(text = element_text(size=17), legend.position = "none")

  
  beta_custom_norm_NMDS(ps.f,color = "Assotiation", shape = "Repeat")

  beta_custom_norm_NMDS_only_bray <- function(ps,  normtype="vst", color="Repeats", shape="Repeats"){
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
    
    ordination.b <- ordinate(ps.varstab, "NMDS", "bray")

    
    #plotting
    p = plot_ordination(ps, ordination.b, type="sample", color, shape, title="NMDS - Bray", 
                        axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
    p <- p + stat_ellipse( type="norm", alpha=0.7)
    
  
    return(p)
  }  

  beta_custom_norm_NMDS_only_bray(ps.f)  
  
  ordination.b <- ordinate(ps.f, "PCoA", "bray")
  p = plot_ordination(ps.f, ordination.b, type="sample", color="Repeats", title="PCoA - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p <- p + stat_ellipse( type="norm", alpha=0.7)
p
  permanova.16cell <- function(ps, dist = "bray"){
    require(phyloseq)
    require(vegan)
    dist <- phyloseq::distance(ps, dist)
    metadata <- as(sample_data(ps), "data.frame")
    ad <- adonis2(dist ~ Repeats, data = metadata, permutations = 10000)
    return(ad)
  }
  
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



