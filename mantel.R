library(phyloseq)
library(vegan)


ps.gr.Arhaea <- subset_taxa(ps.gr, Kingdom=="Archaea")

gr.map.d <- as.data.frame(gr.map.num)
gr.map.num <- mapply(gr.map, FUN=as.numeric)

mantel.res <- bioenv(otu_table(ps.gr.Arhaea), gr.map.d, method="pearson")
summary(mantel.res)

dist.ps.gr <- distance(ps.gr, "bray")
env.dist <- vegdist(scale(gr.map.num[,9]), "euclid")
mantel(dist.ps.gr, env.dist, method="pearson", permutations = 9999)

Job1 = mcparallel(code1())
JobResult1 = mccollect(Job1)


map.red <- map.d[which(map.d$Plant=="s7307"),]

mantel.local <- function(ps, plant, map.d){
    ps.pruned <- prune_samples(sample_data(ps)$Plant %in% c(plant), ps)
    dist <- distance(ps.pruned, "bray")
    map <- map.d[which(map.d$Plant=="s7307"),]
    env.dist <- vegdist(scale(map$pH) , "euclid")
    mant <- mantel(dist, env.dist, method="pearson", permutations = 9999)
    return(mant)
}

map.d %>%
    group_by(Plant,Al) %>%
    summarize(mean_size = mean(BineSeedM, na.rm = TRUE))

map.d %>%
    group_by(Plant,Al,Inoculation) %>%
    summarize(mean_size = cor(SoilAl,BineAl))

cor(map.num$Plant, map.num$BineAl)

veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

otus.ps.vegan <- veganifyOTU(ps)
otus.ps.vegan.descrnames <- otus.ps.vegan
rownames(otus.ps.vegan.descrnames) <- ps@sam_data$Description

vare.cca <- cca(otus.ps.vegan ~ pH + SoilAl, data=map.env)
plot(vare.cca,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,1.5),scaling=2)
points(vare.cca,disp="sites",pch=21,col="red",bg="red",cex=1.3)
text(vare.cca,"sites",pos=3,axis.bp=TRUE)


otus.ps.vegan <- veganifyOTU(ps.7307)
otus.ps.vegan.descrnames <- otus.ps.vegan
rownames(otus.ps.vegan.descrnames) <- ps.7307@sam_data$Description
metadata <- as(sample_data(ps.7307), "data.frame")
vare.cca <- cca(otus.ps.vegan.descrnames ~ pH + SoilAl, data=metadata)
plot(vare.cca,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,1.5),scaling=2)
points(vare.cca,disp="sites",pch=21,col="red",bg="red",cex=1.3)
text(vare.cca,"sites",pos=3,axis.bp=TRUE)

otus.ps.vegan <- veganifyOTU(ps.7307)
otus.ps.vegan.descrnames <- otus.ps.vegan
rownames(otus.ps.vegan.descrnames) <- ps.7307@sam_data$Description
metadata <- as(sample_data(ps.7307), "data.frame")
vare.cca <- cca(otus.ps.vegan.descrnames ~ pH + SoilAl, data=metadata)
plot(vare.cca,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,4),scaling=2)
points(vare.cca,disp="sites",pch=21,col="red",bg="red",cex=1)
text(vare.cca,"sites",pos=4,axis.bp=TRUE)









mantel.local <- function(ps, ){
    gr.map.d <- as.data.frame(gr.map.num)
    gr.map.num <- mapply(gr.map, FUN=as.numeric)
    ps.pruned <- prune_samples(sample_data(ps)$Plant %in% c(plant), ps)
    dist <- distance(ps.pruned, "bray")
    map <- map.d[which(map.d$Plant=="s7307"),]
    env.dist <- vegdist(scale(map$pH) , "euclid")
    mant <- mantel(dist, env.dist, method="pearson", permutations = 9999)
    return(mant)
}