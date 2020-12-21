alignment <- AlignSeqs(DNAStringSet(Al.dna), anchor=NA,verbose=FALSE, processors = 4)
writeXStringSet(alignment, file="align.fasta")

#iqtree -s align.rdp.fasta -nt AUTO -ntmax 40 -bb 1000 -bnni 


# dome transformations of table to Al scikit 

diagdds = phyloseq_to_deseq2(ps.f, ~ Description)                  

diagdds = estimateSizeFactors(diagdds, type="poscounts")
diagdds = estimateDispersions(diagdds, fitType = "local") 
pst <- varianceStabilizingTransformation(diagdds)
pst.dimmed <- t(as.matrix(assay(pst))) 

ps.varstab <- ps.f
otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 

ps.merged.varstab <- merge_samples(ps.varstab, "Description", median)

write.table(ps.merged.varstab@sam_data, file = "metadata_for_rf_merged.txt", sep= "\t", col.names = NA, quote=FALSE)


#for plates dataset

library(phyloseq)
library(tibble)
library(ggplot2)
library(ggrepel)

taxa <- as.matrix(read.csv("taxa.csv" , header=TRUE, sep="\t", row.names = "ID"))


otus <- t(read.csv("intable.otus.csv" , row.names = "ID", sep="\t"))

otus <- t(read.csv("intable.otus.csv" , row.names = "ID", sep="\t"))

mapp <- read.csv("fake_meta.csv" , header=TRUE, sep="\t")
head(mapp)
fake.map <- data.frame(row.names="ID", mapp)

ps <- phyloseq(otu_table((otus), taxa_are_rows=FALSE), 
               sample_data(fake.map), 
               tax_table(taxa))



ordination.b <- ordinate(ps, "NMDS", "bray", type="taxa")
plot_ordination(ps, ordination.b, type="species", title="NMDS - Bray - major OTUs", axes = c(1,2), label = "name", color= "Phylum"  ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 2) 

  
