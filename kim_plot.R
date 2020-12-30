library(vegan)
library(ape)
library(phyloseq)

function(ps, fasta, color){dist <- dist.dna(fasta)
fasta.pcoa <- pcoa(dist, correction="none")
p <- plot_ordination(ps, fasta.pcoa, type = "samples", axes = 1:2,color="kek")
return(p)
}

metadata <- as(sample_data(kim.ps), "data.frame")
adonis2(dist ~ kek, data = metadata)

mapp <- read.csv("map_nodx.csv" , header=TRUE, sep="\t")
kim.map.nodx <- data.frame(row.names="ID", mapp)
kim.nodx.fasta <- read.dna("nodX_art.fas", format = "fasta", as.matrix = TRUE)

kim.fake.otu <- matrix(rexp(200, rate=.1), ncol=54)
colnames(kim.fake.otu)<-rownames(kim.map.nodx)
rownames(kim.fake.otu)<- c("otu1","otu2","otu3","otu4")
kim.ps <- phyloseq(otu_table(t(kim.fake.otu), taxa_are_rows=FALSE), sample_data(kim.map.nodx))

kim_dist_plot <- function(ps,  fasta, color, label){dist <- dist.dna(fasta)
fasta.pcoa <- pcoa(dist, correction="none")
p <- plot_ordination(ps, fasta.pcoa, type = "samples", axes = 1:2,color=color, label=label)
return(p)
}


permanova <- function(ps, Factor, dist = "bray"){
    require(phyloseq)
    require(vegan)
    dist <- distance(ps, dist)
    metadata <- as(sample_data(ps), "data.frame")
    ad <- adonis2(dist ~ Factor, data = metadata)
    return(ad)
}
    
permanova.local <- function(ps, plant, factor){
    ps.pruned <- prune_samples(sample_data(ps)$Plant %in% c(plant), ps)
    dist <- distance(ps.pruned, "bray")
    metadata <- as(sample_data(ps.pruned), "data.frame")
    ad <- adonis2(dist ~ Site, data = metadata)
    return(ad)
}