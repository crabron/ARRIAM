
#splitted up dataset by plant type(genotype)

s1903.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s1903"), ps.f)
s7307.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s7307"), ps.f)
s8353.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s8353"), ps.f)
s8473.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s8473"), ps.f)

#making beta diversity plots

source("storage/scripts_git_done/beta_for_Al.R")
p.1903 <- beta_for_Al(s1903.ps)

permanova <- function(ps, factor = Repeats, distance_metric = "bray"){
    require(vegan)
    require(phyloseq)
    dist <- distance(ps@ , distance_metric)
    metadata <- as(sample_data(ps), "data.frame")
    ad <- adonis2(dist ~ factor, data = metadata)
    return(ad)
}

permanova(s1903.ps, factor=Al, distance_metric="bray")


library(tibble)
library(phyloseq)
otu <- read.csv("OTU.csv" , header=TRUE, sep="\t")
otu <- column_to_rownames(otu, "ID")
otu <- t(as.matrix(otu))


mapp <- read.csv("map_with.csv" , header=TRUE, sep="\t")
map <- column_to_rownames(mapp, "ID")
map <- as.data.frame(t(map))


taxa <- read.csv("taxa.csv" , header=TRUE, sep="\t")
taxa <- column_to_rownames(taxa, "ID")
taxa <- as.matrix(taxa)

ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
             sample_data(map), 
             tax_table(taxa))
             
             

library(phyloseq)             
library(vegan)

# функция преобразует otu-table из phyloseq объекта в otu-table необходимый для пакета vegan
veganifyOTU <- function(physeq){
  require(phyloseq)
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#вместо ps необходимый объект класса phyloseq
otus.ps.vegan <- veganifyOTU(ps)
metadata <- as(sample_data(ps), "data.frame")
#собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
vare.cca <- cca(otus.ps.vegan ~ C_N + Carbon + Ntot + P2O5 + K2O + N_NH4 + N_NO3 + pHh2o + pHCaCl2 + EA + HA , data=metadata)
#рисовалка биплота
plot(vare.cca,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,4),scaling=2)
#стрелочек
points(vare.cca,disp="sites",pch=21,col="red",bg="red",cex=1)
#надписей
text(vare.cca,"sites",pos=4,axis.bp=TRUE)



otus.ps.vegan <- veganifyOTU(ps)
otus.ps.vegan.descrnames <- otus.ps.vegan
rownames(otus.ps.vegan.descrnames) <- ps@sam_data$Description
metadata <- as(sample_data(ps), "data.frame")
vare.cca <- cca(otus.ps.vegan.descrnames ~ C_N + Carbon + Ntot + P2O5 + K2O + N_NH4 + N_NO3 + pHh2o + pHCaCl2 + EA + HA, data=metadata)
plot(vare.cca,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,4),scaling=2)
points(vare.cca,disp="sites",pch=21,col="red",bg="red",cex=1)
text(vare.cca,"sites",pos=4,axis.bp=TRUE)


otus.ps.vegan <- veganifyOTU(ps)
otus.ps.vegan.descrnames <- otus.ps.vegan
rownames(otus.ps.vegan.descrnames) <- ps@sam_data$Description
metadata <- as(sample_data(ps), "data.frame")
vare.cca <- cca(otus.ps.vegan.descrnames ~ C_N + Carbon + Ntot + P2O5 + K2O + N_NH4 + N_NO3 + pHh2o + pHCaCl2 + EA + HA , data=metadata)
plot(vare.cca,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,4),scaling=2)
points(vare.cca,disp="sites",pch=21,col="red",bg="red",cex=1)
text(vare.cca,"sites",pos=4,axis.bp=TRUE)

#вместо ps необходимый объект класса phyloseq
otus.ps.vegan <- veganifyOTU(ps)
metadata <- as(sample_data(ps), "data.frame")
#собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
vare.dbrda <- dbrda(otus.ps.vegan ~ C_N + Carbon + Ntot + P2O5 + K2O + N_NH4 + N_NO3 + pHh2o + pHCaCl2 + EA + HA + Lake_dist, data=metadata, distance = "bray")
#рисовалка биплота
plot(vare.dbrda,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,4),scaling=2)
#стрелочек
points(vare.dbrda,disp="sites",pch=21,col="red",bg="red",cex=1)
#надписей
text(vare.dbrda,"sites",pos=4,axis.bp=TRUE)


vare.cca
vare.cca <- cca(otu ~pHh2o, data=metadata)
mod <- ordistep(little.cca, scope = formula(vare.cca)) little.cca - cca по pHh2o)
anova(mod)
anova(mod, by="terms")
anova(mod, by="mar")

ps.f.dec <-  subset_samples(ps.f.dec, sample_names(ps.f.dec) != 'belimov.64')
ps.f.dec <-  subset_samples(ps.f.dec, sample_names(ps.f.dec) != 'belimov.63')
ps.f.dec <- prune_taxa(taxa_sums(ps.f.dec) > 0, ps.f.dec)

rish <- estimate_richness(ps.f, measures = "Observed")
reads.sum <- as.data.frame(sample_sums(ps.f))
reads.summary <- cbind(rish, reads.sum)

#картиночки  с разными результатами по ридам/otu
p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus.clear, y=log2(reads.clear))) + geom_text(aes(x=otus.clear, y=log2(reads.clear),label=row.names(reads.summary)))
p1 +  geom_point(aes(x=otus.cont, y=log2(reads.cont))) + geom_text(aes(x=otus.cont, y=log2(reads.cont),label=row.names(reads.summary)))

p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus.clear, y=log2(reads.clear))) + geom_point(aes(x=otus.cont, y=log2(reads.cont)),color = "red",alpha=0.5) + geom_segment(aes(x=otus.cont, xend = otus.clear, y = log2(reads.cont), yend = log2(reads.clear)), alpha=0.5,                                                                                                                                                                                                     color = "blue",arrow = arrow(length = unit(0.2, "cm"))) +  geom_text_repel(aes(x=otus.cont,y=log2(reads.cont),label=row.names(reads.summary), size=0.7, alpha=0.5))
p1 +  stat_smooth(aes(x=otus.clear,y=log2(reads.clear)), method = "glm",formula = y ~ ns(x,2))

 ggplot(data=reads.summary) + geom_point(aes(x=otus.clear, y=log2(reads.clear))) + geom_point(aes(x=otus.cont, y=log2(reads.cont)),color = "red",alpha=0.5) + geom_segment(aes(x=otus.cont, xend = otus.clear, y = log2(reads.cont), yend = log2(reads.clear)), alpha=0.5,                                                                                                                                                                                                     color = "blue",arrow = arrow(length = unit(0.2, "cm"))) +  geom_text_repel(aes(x=otus.cont,y=log2(reads.cont),label=row.names(reads.summary), size=2))


ps.f <-  subset_samples(ps.f, sample_names(ps.f) != 'belimov.41')
ps.f <-  subset_samples(ps.f, sample_names(ps.f) != 'belimov.16')
ps.f <- prune_taxa(taxa_sums(ps.f.dec) > 0, ps.f)

s1903.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s1903"), ps.f)
s1903.ps  <- prune_taxa(taxa_sums(s1903.ps ) > 0, s1903.ps )
s7307.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s7307"), ps.f)
s7307.ps <- prune_taxa(taxa_sums(s7307.ps) > 0, s7307.ps)
s8353.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s8353"), ps.f)
s8353.ps <- prune_taxa(taxa_sums(s8353.ps) > 0, s8353.ps)
s8473.ps <- prune_samples(sample_data(ps.f)$Plant %in% c("s8473"), ps.f)
s8473.ps<- prune_taxa(taxa_sums(s8473.ps) > 0,s8473.ps)


otus.ps.vegan <- veganifyOTU(ps.varstab.mod)
metadata <- as(sample_data(ps.varstab.mod), "data.frame")
rownames(otus.ps.vegan) <- metadata$Description
vare.cca <- cca(otus.ps.vegan ~  BineM + NodPerPlant + Mic + pH + NitBine + Nit15 + NitPerPlant + Nit15PerPlant + NitSeed + SoilAl + SoilB + SoilCa + SoilCo + SoilCu + SoilFe + SoilK + SoilMg + SoilMn + SoilMo + SoilNi + SoilP + SoilS + SoilZn , data=metadata)


            
diagdds <- phyloseq_to_deseq2(ps.1903, ~ Description)                  
diagdds = estimateSizeFactors(diagdds, type="poscounts")
diagdds = estimateDispersions(diagdds, fitType = "local") 
pst <- varianceStabilizingTransformation(diagdds)
pst.dimmed <- t(as.matrix(assay(pst))) 
pst.dimmed[pst.dimmed < 0.0] <- 0.0 
ps.varstab.mod <- ps.1903.mod            
otu_table(ps.varstab.mod) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)

otus.ps.vegan <- veganifyOTU(ps.varstab.mod)
metadata <- as(sample_data(ps.varstab.mod), "data.frame")
rownames(otus.ps.vegan) <- metadata$Description
vare.cca <- cca(otus.ps.vegan ~ SoilAl + pH + SoilCo + BineM, data=metadata)

vare.cca
anova(vare.cca, by="terms")
anova(vare.cca, by="mar")

mod <- ordistep(little.cca, scope = formula(vare.cca))
mod$anova



