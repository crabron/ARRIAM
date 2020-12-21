
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


plot_rich_reads <- function(ps.f){
  rish <- estimate_richness(ps.f, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(ps.f))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Description"] <- ps.f@sam_data$Description
  library(ggrepel)
  p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus, y=reads)) + geom_text_repel(aes(x=otus, y=reads,label=rownames(reads.summary))) + theme_bw()
  library(ggforce)
  p1 <- p1  + geom_mark_ellipse(aes(x = otus, y=reads, group = Description, label = Description), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  return(p1)
}
plot_rich_reads(ps.f.dec)
ps.f@sam_data$Description

plot_rich_reads_samlenames <- function(ps.f){
  rish <- estimate_richness(ps.f, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(ps.f))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Repeats"] <- rownames(ps.f@sam_data)
  library(ggrepel)
  p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus, y=reads)) + geom_text_repel(aes(x=otus, y=reads,label=Repeats)) + theme_bw() +
  return(p1)
}
plot_rich_reads_samlenames(ps.f.def)

plot_rich_reads_samlenames_lm <- function(ps.f){
  rish <- estimate_richness(ps.f, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(ps.f))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Description"] <- rownames(ps.f@sam_data)
  reads.summary["Repeats"] <- ps.f@sam_data$Description
  library(ggrepel)
  require(ggforce)
  p1 <- ggplot(data=reads.summary) + geom_point(aes(x=otus, y=reads)) + geom_text_repel(aes(x=otus, y=reads,label=Description)) + theme_bw()+
    geom_smooth(aes(x=otus, y=reads, fill=Repeats),method=lm, se=FALSE) 
    #geom_mark_ellipse(aes(x = otus, y=reads, group = Repeats, label = Repeats, color = Repeats), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  return(p1)
}
plot_rich_reads_samlenames_lm(ps.f.dec)

plot_rich_reads_samlenames_lm <- function(ps.f){
  rish <- estimate_richness(ps.f, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(ps.f))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Description"] <- rownames(ps.f@sam_data)
  reads.summary["Repeats"] <- ps.f@sam_data$Description
  reads.summary["Plant"] <- ps.f@sam_data$Plant
  library(ggrepel)
  require(ggforce)
  p1 <- ggplot(data=reads.summary) + geom_point(aes(y=otus, x=reads, color=Repeats),size=3) + geom_text_repel(aes(y=otus, x=reads, label=Repeats)) + theme_bw()+
    geom_smooth(aes(y=otus, x=reads, fill=Repeats, color=Repeats,),method=lm, se=FALSE, ymin = 1)  
  # geom_mark_ellipse(aes(y = otus, x=reads, group = Repeats, label = Repeats, color = Repeats), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  return(p1)
}
plot_rich_reads_samlenames_lm(ps.f.dec)
p.rich.wc <- plot_rich_reads_samlenames_lm(ps.f)
p.rich.wc
p.rich.dec <- plot_rich_reads_samlenames_lm(ps.f.dec)
p.rich.dec
p.rich.wc + p.rich.dec
# start to making my superior happy
ggarrange(p.rich.wc, p.rich.dec, nrow = 1)

physeq <- ps.f.dec

rish <- estimate_richness(physeq, measures = "Observed")
reads.sum <- as.data.frame(sample_sums(physeq))
reads.summary <- cbind(rish, reads.sum)
colnames(reads.summary) <- c("otus","reads")
reads.summary["Description"] <- rownames(physeq@sam_data)
reads.summary["Repeats"] <- physeq@sam_data$Description
reads.summary["Plant"] <- physeq@sam_data$Plant
library(ggrepel)
require(ggforce)
library(tidyverse)
lm_pvalue <- function(reads.summary, selected, bywhat){
  reads.summary.1903 <-  reads.summary[reads.summary[[bywhat]] %in% c(selected),]
  fit.1903 <- lm(otus ~ reads ,reads.summary.1903)
  fit.1903.sum <- summary(fit.1903)
  fit.1903.sum$coefficients %>% 
    as.data.frame() %>% 
    slice(2) %>%  
    pull(4) -> res 
  return(res)
}

p1 <- ggplot(data=reads.summary) + geom_point(aes(y=otus, x=reads, color=Repeats),size=3) + geom_text_repel(aes(y=otus, x=reads, label=Repeats)) + theme_bw()+
  geom_smooth(aes(y=otus, x=reads, fill=Repeats, color=Repeats,), method=lm, se=FALSE)  
# geom_mark_ellipse(aes(y = otus, x=reads, group = Repeats, label = Repeats, color = Repeats), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))

library(tidyverse)
lm_pvalue_forplot <- function(selected){
  reads.summary.1903 <-  subset(reads.summary, Repeats %in% c(selected))
  fit.1903 <- lm(otus ~ reads ,reads.summary.1903)
  fit.1903.sum <- summary(fit.1903)
  fit.1903.sum$coefficients %>% 
    as.data.frame() %>% 
    slice(2) %>%  
    pull(4) -> res 
  return(res)
}


amp_rarecurve(
  amp.f,
  stepsize = 1000,
  color_by = "Description",
  facet_by = "Plant",
  facet_scales = "fixed"
)

source("~/storage/scripts_git_done/functions.R")
amp.f <- phyloseq_to_amp_without_tree(ps.f)

library(phyloseq)
library(dplyr)
ps.f@sam_data$Inoculation
old_map <- ps.f@sam_data 
levels(old_map) <- c(levels(old_map), "pos")
levels(old_map) <- c(levels(old_map), "neg")
old_map[old_map == "+"] <- "pos"
old_map[old_map == "-"] <- "neg"
old_map

ps.f.ch <- ps.f
mapp <- read.csv("map.tsv" , header=TRUE, sep="\t")
map <- data.frame(row.names="ID", mapp)
sam_data(ps.f.ch) <- map
ps.f.ch@sam_data

setwd("~/storage/al_R/")
write.table(ps.f@sam_data, file = "map.tsv", sep= "\t", col.names = NA, quote=FALSE)

ps.f.s1903 <- prune_samples(sample_data(ps.f.ch)$Plant %in% c("s1903"), ps.f.ch)
ps.f.s1903 <- prune_taxa(taxa_sums(ps.f.s1903) > 0,ps.f.s1903)

ps.f.s7307 <- prune_samples(sample_data(ps.f.ch)$Plant %in% c("s7307"), ps.f.ch)
ps.f.s7307 <- prune_taxa(taxa_sums(ps.f.s7307) > 0,ps.f.s7307)

ps.f.s8473 <- prune_samples(sample_data(ps.f.ch)$Plant %in% c("s8473"), ps.f.ch)
ps.f.s8473 <- prune_taxa(taxa_sums(ps.f.s8473) > 0,ps.f.s8473)

ps.f.s8353 <- prune_samples(sample_data(ps.f.ch)$Plant %in% c("s8353"), ps.f.ch)
ps.f.s8353 <- prune_taxa(taxa_sums(ps.f.s8353) > 0,ps.f.s8353)

physeq <- ps.f.s1903
physeq@sam_data$Al

physeq.Al.pos <- prune_samples(sample_data(physeq)$Al %in% c("pos"), physeq)
physeq.Al.pos <- prune_taxa(taxa_sums(physeq.Al.pos) > 0, physeq.Al.pos)

physeq.Al.neg <- prune_samples(sample_data(physeq)$Al %in% c("neg"),physeq)
physeq.Al.neg <- prune_taxa(taxa_sums(physeq.Al.neg) > 0, physeq.Al.neg)

physeq.In.pos <- prune_samples(sample_data(physeq)$Inoculation %in% c("pos"), physeq)
physeq.In.pos <- prune_taxa(taxa_sums(physeq.In.pos) > 0,physeq.In.pos)

physeq.In.neg <- prune_samples(sample_data(physeq)$Inoculation %in% c("pos"), physeq)
physeq.In.neg <- prune_taxa(taxa_sums(physeq.In.neg) > 0,physeq.In.neg)

physeq@sam_data


beta_custom_norm_NMDS_elli <- function(ps, seed = 7888, normtype="vst", color="Repeats"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  library(ggforce)
  
  
  # beta_NMDS <- function(){
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Description)                  
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
  p1 = plot_ordination(ps, ordination.b, type="sample", color="Description", title="NMDS - Bray", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Description, label = Description), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  p2 = plot_ordination(ps, ordination.u, type="sample", color="Description", title="NMDS - unifrac", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Description, label = Description), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  p3 = plot_ordination(ps, ordination.w, type="sample", color="Description", title="NMDS - wunifrac", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) + 
    geom_mark_ellipse(aes(group = Description, label = Description), label.fontsize = 10, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1, p2, p3, ncol = 3 , nrow = 1, common.legend = TRUE, legend = "none", font.label = list(size = 12, face = "bold", color ="black"))
  
  return(p.all)
}

p.beta <- beta_custom_norm_NMDS_elli(physeq)

alpha.custom.Al <- function(ps.f.i, arrga = "PD"){
  require(picante)
  pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)

  pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
  pd <- cbind(pd, ps.f.i@sam_data, pd1 )

  bmi <- levels(as.factor(pd$Al))
  bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
  p1 <- ggviolin(pd, x = "Al", y = arrga,
                 add = "boxplot", fill = "Al") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
    scale_x_discrete(limits = c("pos","neg")) + theme(axis.title.x = element_blank(), legend.title = element_blank())
  return(p1)
}


alpha.all.Al <- function(physeq.max, physeq.min , alpha.custom, arrga = "PD")  {

  p.max.o <- alpha.custom(physeq.max, arrga = "Observed")
  p.max.sh <- alpha.custom(physeq.max, arrga = "Shannon")
  p.max.is <- alpha.custom(physeq.max, arrga = "InvSimpson")
  p.max.pd <- alpha.custom(physeq.max, arrga = "PD")

  p.min.o <- alpha.custom(physeq.min , arrga = "Observed")
  p.min.sh <- alpha.custom(physeq.min , arrga = "Shannon")
  p.min.is <- alpha.custom(physeq.min , arrga = "InvSimpson")
  p.min.pd <- alpha.custom(physeq.min , arrga = "PD")

  p.alpha.Al <- ggarrange(p.max.o, p.max.sh, p.max.is,  p.max.pd, p.min.o, p.min.sh, p.min.is, p.min.pd, ncol = 4 ,label.x = 0.105, nrow = 2, common.legend = TRUE)
  return(p.alpha.Al)
}
alpha.custom.Al(physeq.In.pos, arrga = "Observed")
p.alpha.Al <- alpha.all.Al(physeq.In.pos, physeq.In.neg, alpha.custom.Al)

p.alpha.Al

Des.soil.w.simper.Al <- function(ps){
  require(DESeq2)
  require(vegan)
  require(tibble)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Al)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$Al), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
} 

des.Al <- physeq.Al
physeq.Al <-  Des.soil.w.simper.Al(physeq)

filter.des <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <- des[(des$log2FoldChange > 0),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}


des.Al.filt <- filter.des(des.Al, alpha = 0.1, baseMean = 10)
length(rownames(des.Al.filt))
View(des.Al)

Des.soil.w.simper.In <- function(ps){
  require(DESeq2)
  require(vegan)
  require(tibble)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Al)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$In), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
} 

des.In  <-  Des.soil.w.simper.In(physeq)

filter.des <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <- des[(des$log2FoldChange > 0),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}
des.In.filt <- filter.des(des.In, alpha = 0.1, baseMean = 10)
length(rownames(des.In.filt))
View(des.In)
# for pod soil



des.black.dnaextc <-  Des.soil.w.simper(ps.s.nosilt.pod.dnaextc)
des.pod.ex.filtered <- filter.des.out.black.dnaextc(des.black.dnaextc, alpha = 0.1, baseMean = 50)
length(rownames(des.black.ex.filtered))
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
  diagdds = phyloseq_to_deseq2(ps, ~ Description) # вставить своё                 
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


# WGCNA
library(vegan)
library(WGCNA)
physeq <- ps.f.dec
sample_names(physeq)
sample_names(physeq) <- paste0(ps.f.dec@sam_data$Description, "_", stringr::str_split_fixed(sample_names(physeq), "elimov.", 2)[,2])
sample_names(physeq) %in% c("s8473_In_41", "s8473_In_Al_45")
physeq <- prune_samples(sample_names(physeq) != c("s8473_In_41", "s8473_In_Al_45"), physeq)
physeq <- prune_taxa(taxa_sums(physeq) > 10, physeq)
physeq@sam_data@row.names
physeq
physeq.var <- varstab(physeq)
otus <- veganifyOTU(physeq.var)
head(otus)
gsg = goodSamplesGenes(otus.d[-1], verbose = 3)
gsg
otus.d <- as.data.frame(otus)

sampleTree = hclust(dist(otus.d), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
metadata <- physeq@sam_data

OTUSamples = rownames(otus.d)
MetaSamples = rownames(metadata)
traitRows = match(OTUSamples, MetaSamples)
traitRows
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage()

sampleTree2 = hclust(dist(otus.d), method = "average")
traitColors = numbers2colors(metadata[,5:111], signed = FALSE)
length(metadata)
View(data.frame(metadata))
metadata.little <- metadata
# small metadata?
enableWGCNAThreads()
powers = c(c(1:10), seq(from = 11, to=30, by=1))
sft = pickSoftThreshold(otus.d, powerVector = powers, verbose = 5, networkType = "signed")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.8,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower = 17;
adjacency = adjacency(otus.d, power = softPower, type = "signed");
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM
TaxaTree = hclust(as.dist(dissTOM), method = "average");


TaxaTree = hclust(as.dist(dissTOM), method = "average");
  sizeGrWindow(12,9)
plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
minModuleSize = 50;
dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Taxa dendrogram and module colors")

MEList = moduleEigengenes(otus.d, colors = dynamicColors)

MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.30
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(otus.d, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
mergedColors
plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
save(MEs, moduleLabels, moduleColors, TaxaTree, file = "Monterey-networkConstruction-stepByStep.RData")
nTaxa = ncol(otus.d);
nSamples = nrow(otus.d);
MEs0 = moduleEigengenes(otus.d, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, metadata, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metadata),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

colnames(met.d)
length(colnames(met.d))
met.d <- data.frame(metadata)
met.d.num <- metadata[,6:110]
met.d.num.scaled <- scale(met.d.num)
View(met.d.num.scaled)
# nmds for factors?
dist.met.scaled <- vegdist(t(met.d.num.scaled), method = "euclidean", na.rm = TRUE)
dist.met.scaled
met.nmds <-
  metaMDS(t(met.d.num.scaled),
          distance = "euclidean",na.rm = TRUE,
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
dev.off()
stressplot(met.nmds)
plot(met.nmds, "sites")
orditorp(met.nmds, "sites")
# scaled - shit
# unscaled - not work at all

res <- pcoa(dist.met.scaled)
res$values
biplot(res)

# pcoa -same shit
hclust.met <-  hclust(dist.met.scaled)
plot(hclust.met, main = "metadata scaled euclidian")

# 
sampleTree2 <-  hclust(dist(met.d.num), method = "average")
traitColors = numbers2colors(met.d.num, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(met.d.num),
                    main = "Sample dendrogram and trait heatmap")
View(met.d)
physeq <- physeq.var
permanova.forloop <- function(factor){
  require(vegan)
  require(phyloseq)
  dist <- phyloseq::distance(physeq, "bray")
  metadata <- as(sample_data(physeq@sam_data), "data.frame")
  ad <- adonis2(as.formula(paste( "dist ~", "Plant/", factor)), data = metadata)
  ad <- cbind2(ad[2,][3], ad[2,][5])
  return(ad)
}

plot_rich_reads_samlenames_lm(physeq)


dist <- phyloseq::distance(physeq, "bray")
metadata <- as(sample_data(physeq@sam_data), "data.frame")
ad <- adonis2(dist ~ pH, data = metadata)
View(ad)
class(ad)
cbind2(ad[1,][3], ad[1,][5])
colnames(physeq@sam_data)[7:8]
permanova.forloop("pH")
lapply(colnames(physeq@sam_data), permanova.forloop)
length(metadata)
metadata <- as(sample_data(physeq@sam_data), "data.frame")
metadata
colnames(physeq@sam_data)[6:110]
lapply(colnames(physeq@sam_data)[4:8], permanova.forloop)
colnames(physeq@sam_data)
?as.formula

saveRDS(ps.f.dec, file = "~/storage/al_R/rf/ps.f.dec.rds")
readRDS(file = "~/storage/al_R/rf/ps.f.dec.rds")
