dada2_pipe <- function(path = "in/", path_trein_set = "data/silva_nr_v132_train_set.fa", path_trein_set_species = "data/silva_species_assignment_v132.fa", map_path, name_Brief, truncLen = "220,180", maxEE = "2,5", mult = TRUE, mlt = NULL){
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
      
      out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180), trimLeft=20, maxN=0, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, multithread=mult)
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
      
      system("conda activate qiime2-2019.1
      sed -i '1s/^/#OTU_ID/' otu_table.txt 
      biom convert -i otu_table.txt -o otu_table.biom --to-hdf5
      qiime tools import \
      --input-path rep_seq.fasta \
      --output-path rep_seq.qza \
      --type 'FeatureData[Sequence]'
      qiime fragment-insertion sepp --i-representative-sequences rep_seq.qza  --o-tree insertion-tree.qza  --o-placements insertion-placements.qza --p-threads 20
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
      
      mapp <- read.csv(map_path , header=TRUE, sep="\t")
      map <- data.frame(row.names="ID", mapp)
      
      ps <- phyloseq(otu_table(filt.otu.matrix, taxa_are_rows=FALSE), 
                     sample_data(map), 
                     tax_table(taxa),
                     phy_tree(tree))
      return(ps)
}

pop.taxa <- function(physeq, badTaxa){
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



