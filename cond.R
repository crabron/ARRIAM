# Svir dataset
# 5-7 sites
# The first approaches -- ls(SS) workshop

func.more <- function(x){ 
  ress <- x[1] > x[2] | x[1] > x[3] | x[3] > x[4]
  return(ress)
}
otus.lapp <- as.data.frame(apply(otus.t@.Data, 1, func.more))

colnames(otus.lapp) <- c("Cond")
View(otus.lapp)
View(otus.t@.Data)
head(otus.t)

library(phyloseq)             
library(vegan)

# функция преобразует otu-table из phyloseq объекта в otu-table необходимый для пакета vegan
veganifyOTU <- function(physeq){
  require(phyloseq)
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#вместо ps необходимый объект класса phyloseq
otus.ps.vegan <- veganifyOTU(ps.s.merged)
metadata <- as(sample_data(ps.s.merged), "data.frame")
#собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
vare.cca <- cca(otus.ps.vegan ~ Mean_depth + N.total + Р2О5 + К2О , data=metadata)

colnames(merged.metadata)

svir.ggcca.species <- function(vare.cca){
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  require(ggvegan)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.8,
                                                                                                                                                                                                             color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
}


ps.s.5 <- prune_samples(sample_data(ps.s)$Site %in% c("5"), ps.s)
ps.s.bf <- prune_samples(sample_data(ps.s)$Horizont %in% c("BF"), ps.s)
ps.s.merged.5 <- prune_samples(sample_data(ps.s.merged)$Run %in% c("5"), ps.s.merged)
sample_data(ps.s.merged)

Des.ls <- function(ps){
  require(DESeq2)
  require(phyloseq)
  diagdds = phyloseq_to_deseq2(ps, ~  Site)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Site), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
}  

kate.ggcca.sites <- function(vare.cca){
  require(ggvegan)
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("sites","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "sites"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.8,
                                                                                                                                                                                                         color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
}

kate.ggcca.sites(ps.s)


ps.s.5.pr <- prune_taxa(taxa_sums(ps.s.5) > 10, ps.s.5)

amp_heatmap(amp, group_by = "Site", facet_by = "Horizont", tax_show = 40,tax_aggregate = "Genus", tax_add = "Family", color_vector = c("white","azure4"))



beta_NMDS <- function(ps){
  require(phyloseq)
  require(ggplot2)
  require(DESeq2)

  diagdds <- phyloseq_to_deseq2(ps, ~ Horizont)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0 
  ps.varstab <- ps            
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)
  ps	<- ps.varstab
  
  ordination.b <- ordinate(ps, "NMDS", "bray")
  p = plot_ordination(ps, ordination.b, type="sample", color="Horizont", shape = "Site", title="NMDS - bray", axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  return(p1.1)
}

ps.s.5.oe.pr <- prune_samples(sample_data(ps.s.5.pr)$Horizont %in% c("O", "E"), ps.s.5.pr)
des.horisont.nh4 <- Des.ls(ps.s.5.oe.pr)
length(which(des.horisont$padj > 0.05))
length(which(des.horisont$padj < 0.05))
length(which(des.horisont.nh4$padj > 0.05))
length(which(des.horisont.nh4$padj < 0.05))

ps.s.5.pr <- prune_taxa(taxa_sums(ps.s.5) > 100, ps.s.5)

library(DESeq2)

diagdds <- phyloseq_to_deseq2(ps.s.5.pr, ~ Horizont )                  
diagdds = estimateSizeFactors(diagdds, type="poscounts")
diagdds = estimateDispersions(diagdds, fitType = "local") 
pst <- varianceStabilizingTransformation(diagdds)
ps.varstab <- t(as.matrix(assay(pst))) 
ps.varstab <- ps.s.5.pr
otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 

cbind.tables <- function(des.horisont){
  temp.df.more <- des.horisont[des.horisont$padj > 0.05,]
  table.df.more <- table(temp.df.more$Phylum)
  sum.more <- sum(table.df.more)
  table.df.more <- table.df.more/sum.more
 temp.df.less <- des.horisont[des.horisont$padj < 0.05,]
  table.df.less <- table(temp.df.less$Phylum)
  sum.less <- sum(table.df.less)
  table.df.less <- table.df.less/sum.less
  table.df <- as.data.frame(cbind(table.df.less, table.df.more))
  table.df <- table.df[rowSums(table.df) > 0,] 
  return(table.df)
}

View(table.df)
table.df <- cbind.tables(des.oe)

library(tibble)

table.fam.all <- rownames_to_column(.data = table.df, var = "Phylum")
dfm <- melt(table.df, id=c(table.fam.all$table.df.less, table.fam.all$table.df.more))
dfm <- cbind(dfm, table.fam.all$Phylum)
p <- ggplot(data=dfm, aes(x=dfm$`table.fam.all$Phylum`, y=value)) + geom_col(aes(fill = variable), position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- ggplot(dfm,aes(x = dfm$`table.fam.all$Family`,y = value)) +   geom_col(aes(fill = variable),stat = "identity",position = "dodge") 

#make DeSeq2 table from bf for 5-6-7
#ps.s.bf.pr.100 <- prune_taxa(taxa_sums(ps.s.bf) > 100, ps.s.bf)
Des.tr <- function(ps.s.bf){
  des.bf.all <- Des.ls(ps.s.bf)
  des.bf.all.2 <- des.bf.all
  des.bf.all.2 <- des.bf.all.2[des.bf.all.2$padj < 0.05,] 
  des.bf.all.2 <- des.bf.all.2[des.bf.all.2$baseMean > 30,]
  des.bf.all.2 <- des.bf.all.2[!is.na(des.bf.all.2$padj),]
  return(des.bf.all.2)
}

Des.ch <- function(ps){
  require(qgraph)
  des.all <- Des.tr(ps)
  ps.meangfl <- prune_taxa(rownames(des.all), ps)
  ps.5.6 <- prune_samples(sample_data(ps.meangfl)$Site %in% c("5", "6"), ps.meangfl)
  ps.6.7 <- prune_samples(sample_data(ps.meangfl)$Site %in% c("6", "7"), ps.meangfl)
  des.5.6 <- Des.ls(ps.5.6)
  des.6.7 <- Des.ls(ps.6.7)
  des.res <- cbind(des.5.6$log2FoldChange, des.6.7$log2FoldChange)
  rownames(des.res) <- rownames(des.all)
  des.res[is.na(des.res) == TRUE] <- 0
  Q <- qgraph(cor(t(des.res)), layout = "spring")
  return(des.res)
}



filter.des.tables <- function(des.horisont){
  temp.df.more <- des.horisont[des.horisont$padj > 0.05,]
  table.df.more <- table(temp.df.more$Phylum)
  sum.more <- sum(table.df.more)
  table.df.more <- table.df.more/sum.more
  temp.df.less <- des.horisont[des.horisont$padj < 0.05,]
  table.df.less <- table(temp.df.less$Phylum)
  sum.less <- sum(table.df.less)
  table.df.less <- table.df.less/sum.less
  table.df <- as.data.frame(cbind(table.df.less, table.df.more))
  table.df <- table.df[rowSums(table.df) > 0,] 
  return(table.df)
}


require(DESeq2)
require(phyloseq)
diagdds = phyloseq_to_deseq2(ps.s.bf, ~  Site)                  
diagdds = estimateSizeFactors(diagdds, type="poscounts")
diagdds = estimateDispersions(diagdds, fitType = "local") 
diagdds = DESeq(diagdds)

write.table(des.bf.all, file = "~/storage/svir/LS/not.nice.reads/some.results/des.bf.all.tsv", sep= "\t", col.names = NA, quote=FALSE)


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
  p1 = plot_ordination(ps, ordination.b,color="Horizont", shape="Site", type="sample", title="NMDS - Bray", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Repeats, label = Repeats), label.fontsize = 12, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  p2 = plot_ordination(ps, ordination.u,color="Horizont", shape="Site", type="sample", title="NMDS - unifrac", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) +
    geom_mark_ellipse(aes(group = Repeats, label = Repeats), label.fontsize = 12, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  p3 = plot_ordination(ps, ordination.w,color="Horizont", shape="Site", type="sample", title="NMDS - wunifrac", 
                       axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) + 
    geom_mark_ellipse(aes(group = Repeats, label = Repeats), label.fontsize = 12, label.buffer = unit(2, "mm"), label.minwidth = unit(5, "mm"),con.cap = unit(0.1, "mm"))
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1, p2, p3, ncol = 3 , nrow = 1, common.legend = TRUE, legend = "none", font.label = list(size = 12, face = "bold", color ="black"))
  
  return(p.all)
}
library(ggpubr)
p1 <- beta_custom_norm_NMDS_elli(ps.s)
p2 <- beta_custom_norm_NMDS_elli(ps.f.chr)
ggarrange(p1, p2, ncol = 1, nrow = 2, label.y = 0.99, label.x = 0.8,  labels = c("Дерново-подзолистая почва","Чернозем"), font.label = list(size = 14, face = "bold", color ="black"))

ps.s.bf <- prune_samples(sample_data(ps.s)$Horizont %in% c("BF"), ps.s)
ps.s.bf <- prune_taxa(taxa_sums(ps.s.bf) > 0, ps.s.bf) 

ps.s.bf.57 <- prune_samples(sample_data(ps.s.bf)$Site %in% c("5","7"), ps.s.bf)
ps.s.bf.57 <- prune_taxa(taxa_sums(ps.s.bf.57) > 0, ps.s.bf.57) 

ps.s.c <- prune_samples(sample_data(ps.s)$Horizont %in% c("C"), ps.s)
ps.s.c <- prune_taxa(taxa_sums(ps.s.c) > 0, ps.s.c) 

ps.s.c.57 <- prune_samples(sample_data(ps.s.c)$Site %in% c("5","7"), ps.s.c)
ps.s.c.57 <- prune_taxa(taxa_sums(ps.s.c.57) > 0, ps.s.c.57) 

tree <- read_tree(treefile="tree_s.nwk")


ps <- phy_tree(tree)
ps
phy_tree(ps.s) <- tree
ps.s

,color="Horizont", shape="Site",

Des.soil.w.simper <- function(ps){
  require(DESeq2)
  require(vegan)
  require(tibble)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Site)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$Site), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
}      

des.bf.sites <-  Des.soil.w.simper(ps.s.c.57)

View(des.bf.sites)


filter.des.out <- function(des, alpha=0.1, baseMean = 15){
  des <- des[(des$padj < alpha), ]
  des <- des[(des$baseMean > baseMean),]
  des <-  des[(is.na(des$padj) != TRUE),]
  return(des)
}

des.otus.n1.pr <- des.otus.n1[rownames(des.otus.chr.n1[(is.na(des.otus.chr.n1$padj) != TRUE),]),]
des.otus.n1.pr <- des.otus.n1.pr[rownames(des.otus.pod.n1[(is.na(des.otus.pod.n1$padj) != TRUE),]),]

des.otus.n1.soils.pr <- filter.des.out(des.bf.sites, alpha = 0.05, baseMean = 75)
length(rownames(des.otus.n1.soils.pr))

ps.f.pr.soils <- prune_taxa(rownames(des.otus.n1.soils.pr), ps.s) 

# ps.f.merged.soil <- merge_samples(ps.f.pr.soils, "Soil")
# melted.soils <- psmelt(ps.f.merged.soil)
# adundance <- melted.soils$Abundance / sum(melted.soils$Abundance)*100
# sigtab = res[(res$padj < alpha), ]
phyloseq <-  ps.f.pr.soils
require(phyloseq)
require(ggtree)
library(scales)
library(ggplot2)
library(ggstance)
library(ape)
library(dplyr)
tree <- phyloseq@phy_tree
taxa.pruned <- as.data.frame(phyloseq@tax_table@.Data)
taxa.pruned <- taxa.pruned %>%  mutate_all(as.character)
taxa.pruned$number <- seq.int(from = nrow(taxa.pruned), to = 1)
# taxa.pruned$number <- seq.int(to = nrow(taxa.pruned), from = 1)
taxa.pruned$taxa <- ifelse(is.na(taxa.pruned$Genus), taxa.pruned$Family, taxa.pruned$Genus)
taxa.pruned$taxa <- ifelse(is.na(taxa.pruned$Family), taxa.pruned$Order, taxa.pruned$Genus)
taxa.pruned[taxa.pruned == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
taxa.pruned$taxa2 <- ifelse(is.na(taxa.pruned$Species), with(taxa.pruned, paste0(taxa)), with(taxa.pruned, paste0(taxa, " ", Species )))
taxa.pruned$taxa3 <- ifelse(taxa.pruned$Phylum == "Proteobacteria", with(taxa.pruned, paste0(taxa.pruned$number, ".", taxa2, " // ", Class)), with(taxa.pruned, paste0(taxa.pruned$number, ".", taxa2, " // ", Phylum)))
class(taxa.pruned$Kingdom)
tree$tip.label <- taxa.pruned$taxa3
#tree$tip.label
# tree <- root(tree, which(tree$tip.label %in% c("49.NA // Crenarchaeota")))
p <- ggtree(tree, ladderize = F) + geom_tiplab(mapping = aes(), align=TRUE, linesize=.5) + xlim(NA, 4)
p

df <- data.frame(id=tree$tip.label, value = des.otus.n1.soils.pr$log2FoldChange)
df$value <- df$value/5
p2 <- facet_plot(p, panel = " Site 5 / Site 7(waterlogged)", data = df, geom = geom_barh, 
                 mapping = aes(value), fill = "Plum", 
                 stat='identity', width = 0.8) 


# p2 <- facet_plot(p2, panel = "Abundance", data = df, geom = geom_barh, 
# mapping = aes(value),
# stat='identity', width = 0.8) 
p2
library(grid)
gt = ggplot_gtable(ggplot_build(p2))
gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.6*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
grid.draw(gt) # plot with grid draw


ps.s <- prune_taxa(taxa_sums(ps.s) > 0, ps.s) 
ps.s
ps.s@sam_data
taxa_s <- as.data.frame(ps.s@tax_table@.Data)
brief <- read.csv("~/storage/svir/new/dada2_pipeline/svir.Brief.tsv" , header=TRUE, sep="\t")
brief.f <- brief %>% rownames_to_column("ID") %>% filter(row.names(brief) %in% row.names(taxa_s))
head(rownames(taxa_s))
length(rownames(brief.f))

library(seqinr)
View(brief.f)
briefToSeq.ls <- as.list(brief.f$x)
briefToSeq.names <- as.list(brief.f$ID)

write.fasta( briefToSeq.ls, briefToSeq.names,  file = "~/storage/svir/new/dada2_pipeline/rep_s.fasta", as.string = FALSE, open = "w")
write.table(t(ps.s@otu_table), file = "~/storage/svir/new/dada2_pipeline/otu_table_s.txt", sep= "\t", col.names = NA, quote=FALSE)
?write.fasta
"
conda activate qiime2-2019.1
sed -i '1s/^/#OTU_ID/' otu_table_s.txt 
biom convert -i otu_table_s.txt -o otu_table_s.biom --to-hdf5
qiime tools import \
--input-path rep_s.fasta \
--output-path rep_s.qza \
--type 'FeatureData[Sequence]'
qiime fragment-insertion sepp --i-representative-sequences rep_s.qza  --o-tree data/insertion-tree.qza  --o-placements data/insertion-placements.qza --p-threads 40
qiime tools import \
--input-path otu_table_s.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path feature-table.qza
qiime fragment-insertion filter-features \
--i-table feature-table.qza \
--i-tree data/insertion-tree.qza \
--o-filtered-table filtered_table.qza \
--o-removed-table removed_table.qza
unzip -p filtered_table.qza */data/* > filtered_table.biom
biom convert -i  filtered_table.biom -o filtered_table_s.txt --to-tsv
sed -i '1d' filtered_table.txt
unzip -p data/insertion-tree.qza */data/* > tree_s.nwk

"
source("~/storage/scripts_git_done/functions.R")


#alpha plots
ps.s.o <- prune_samples(sample_data(ps.s)$Horizont %in% c("O"), ps.s)
ps.s.o <- prune_taxa(taxa_sums(ps.s.o) > 0, ps.s.o) 

ps.s.e <- prune_samples(sample_data(ps.s)$Horizont %in% c("E"), ps.s)
ps.s.e <- prune_taxa(taxa_sums(ps.s.e) > 0, ps.s.e) 

ps.s.bf <- prune_samples(sample_data(ps.s)$Horizont %in% c("BF"), ps.s)
ps.s.bf <- prune_taxa(taxa_sums(ps.s.bf) > 0, ps.s.bf) 



alpha.custom <- function(ps.f.i, arrga = "PD"){
  require(picante)
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)
  
  pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
  pd <- cbind(pd,ps.f.i@sam_data, pd1 )
  
  bmi <- levels(pd$Site)
  bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
  p1 <- ggviolin(pd, x = "Site", y = arrga,
                 add = "boxplot", fill = "Site") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
     theme(axis.title.x = element_blank(), legend.title = element_blank())
  return(p1)
}

p.alpha.Cdt.oo <- alpha.custom(ps.s.e, arrga = "Observed")
p.alpha.Cdt.sh <- alpha.custom(ps.s.e, arrga = "Shannon")
p.alpha.Cdt.is <- alpha.custom(ps.s.e, arrga = "InvSimpson")
p.alpha.Cdt.pd <- alpha.custom(ps.s.e, arrga = "PD")

p.alpha.wt.oo <- alpha.custom(ps.s.o, arrga = "Observed")
p.alpha.wt.sh <- alpha.custom(ps.s.o, arrga = "Shannon")
p.alpha.wt.is <- alpha.custom(ps.s.o, arrga = "InvSimpson")
p.alpha.wt.pd <- alpha.custom(ps.s.o, arrga = "PD")

p.alpha.4.oo <- alpha.custom(ps.s.bf, arrga = "Observed")
p.alpha.4.sh <- alpha.custom(ps.s.bf, arrga = "Shannon")
p.alpha.4.is <- alpha.custom(ps.s.bf, arrga = "InvSimpson")
p.alpha.4.pd <- alpha.custom(ps.s.bf, arrga = "PD")

p.alpha.5.oo <- alpha.custom(ps.s.c, arrga = "Observed")
p.alpha.5.sh <- alpha.custom(ps.s.c, arrga = "Shannon")
p.alpha.5.is <- alpha.custom(ps.s.c, arrga = "InvSimpson")
p.alpha.5.pd <- alpha.custom(ps.s.c, arrga = "PD")

p1 <- ggarrange(p.alpha.Cdt.oo, p.alpha.Cdt.sh, p.alpha.Cdt.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p2 <- ggarrange(p.alpha.wt.oo, p.alpha.wt.sh,p.alpha.wt.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p3 <- ggarrange(p.alpha.4.oo, p.alpha.4.sh,p.alpha.4.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p4 <- ggarrange(p.alpha.5.oo, p.alpha.5.sh,p.alpha.5.pd , ncol = 3 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
ggarrange(p1, p2, p3, p4, ncol = 1 , nrow = 4, label.y = 0.9855, label.x = 0.29,  labels = c("horizon O","horizon E","horizon BF", "horizon C"), font.label = list(size = 14, face = "bold", color ="black"))  

ps.f.i <- ps.s.e
require(picante)
require(phyloseq)
require(ggpubr)
library(ggplot2)

pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)
pd
pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
pd <- cbind(pd,ps.f.i@sam_data, pd1 )
pd
bmi <- levels(pd$Repeats)
bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
bmi.pairs
p1 <- ggviolin(pd, x = "Repeats", y = "Observed",
               add = "boxplot", fill = "Repeats") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
  scale_x_discrete(limits = c( "Site 5","Site 6", "Site 7")) + theme(axis.title.x = element_blank(), legend.title = element_blank())
p1
p1 <- ggviolin(pd, x = "Repeats", y = "Observed", add = "boxplot", fill = "Repeats")
p2 <- p1 + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test")
p2 + theme(axis.title.x = element_blank(), legend.title = element_blank())

ps.s.e
return(p1)


#trying corncob tutorial 

library(corncob)
library(magrittr)
library(phyloseq)
ps.s.bf.57@sam_data
soil <- ps.s.bf.57 %>% tax_glom("Phylum")
soil
tax_table(soil)
corncob <- bbdml(formula = Seq1 ~ 1, phi.formula =  ~ 1,data = soil)
corncob_da <- bbdml(formula = Seq1 ~ Site, phi.formula =  ~ Site,data = soil)
lrtest(corncob, corncob_da)
plot(corncob)
plot(corncob_da, AA=TRUE, color = "Site")
View(corncob)
summary(corncob_da)
fullAnalysis <-  differentialTest(test = "Wald", boot = FALSE,  formula = ~ Site, phi.formula = ~ Site, phi.formula_null = ~ 1, data =soil, formula_null = ~ 1, fdr_cutoff = 0.05)
corn.all <-  differentialTest(test = "Wald", boot = FALSE,  formula = ~ Site, phi.formula = ~ Site, phi.formula_null = ~ 1, data =ps.s.bf.57, formula_null = ~ 1, fdr_cutoff = 0.05)
View(cbind(des.57, corn.all$p_fdr))
head(rownames(corn.all$p_fdr))
length(corn.all$p_fdr)
length(rownames(des.57))
des.filt <- des.bf.sites[corn.all$p_fdr, ]
corn.filt <- as.data.frame(corn.all$p_fdr)[rownames(des.filt),]

library(tidyverse)
library(naturalsort)
des.57 %>% as_tibble(rownames = NA) %>% rownames_to_column() %>% arrange(rowname, naturalsort(rowname)) 
des.wtf <- des.57 %>% as_tibble(rownames = NA) %>% rownames_to_column() %>% 
  dplyr::mutate(order = naturalsort(rowname)) %>% 
  mutate(rowname = factor(rowname, levels = order)) %>% arrange(rowname)

*View(des.wtf)
?arrange
View(corn.all$all_models[1])
View
View(cbind(des.57, corn.all$p_fdr))

#some test with bray-curtis and other beta

View(ps.s.bf.57@tax_table@.Data)
ps.s.bf.57.Sideroxydans <-subset_taxa(ps.s.bf.57, Genus == "Sideroxydans")
ps.s.bf.57.Sideroxydans@otu_table
dist.bf.57.Sideroxydans  <- phyloseq::distance(ps.s.bf.57.Sideroxydans, method="bray", type="taxa")
dist.bf.57.Sideroxydans

bray.steppo <- function(some_number){
  some_ps <- prune_taxa(taxa_sums(some_ps) > some_number , some_ps) 
  otus.dt <- as.data.frame(some_ps@otu_table@.Data)
  jacco.mean <- mean(vegdist(otus.dt, method = "bray"))
  return(jacco.mean)
}

wuni.steppo <- function(some_number){
  some_ps <- prune_taxa(taxa_sums(some_ps) > some_number , some_ps) 
  jacco.mean <- mean(phyloseq::distance(some_ps, method = "wunifrac", type = "samples"))
  return(jacco.mean)
}

uni.steppo <- function(some_number){
  some_ps <- prune_taxa(taxa_sums(some_ps) > some_number , some_ps) 
  jacco.mean <- mean(phyloseq::distance(some_ps, method = "unifrac", type = "samples"))
  return(jacco.mean)
}

some_ps <- ps.s.bf.57
while(which(length(taxa_sums(some_ps)) <= 30)))
seq.list <- seq(min(taxa_sums(some_ps)), max(taxa_sums(some_ps)), by=10)
seq.list <- seq(min(taxa_sums(some_ps)), 500, by=10)
bray.list <- sapply(seq.list, bray.steppo)
wuni.list <- sapply(seq.list, wuni.steppo)
uni.list <- sapply(seq.list, uni.steppo)
plot(bray.list)
plot(wuni.list)
plot(uni.list)

require(dplyr)
require(ggpubr)

ggarrange(p1, p2б, labels = c("bray","wunifrac")) 

# the last version of beta-steppo code

library(tidyverse)

# the anonimeus function - because I can
len.col <- function(x) length(taxa_names(prune_taxa(taxa_sums(some_ps) >= x, some_ps))) >= 30
# ctreate a sequence list to literate throught
seq.list <- seq(0, max(taxa_sums(some_ps)), by=10) %>% purrr::head_while(len.col)
# bray: prune by the phyloseq package functions - veganisation - bray by vegan, simple mean 
bray.steppo <- function(some_number, distance = "bray", ps = some.ps){
  some_ps <- prune_taxa(taxa_sums(some_ps) > some_number , some_ps) 
  otus.dt <- as.data.frame(some_ps@otu_table@.Data)
  jacco.mean <- mean(vegdist(otus.dt, method = distance))
  return(jacco.mean)
}

bray.steppo.nm <- function(some_number, distance = "bray", ps = some.ps){
  some_ps <- prune_taxa(taxa_sums(some_ps) > some_number , some_ps) 
  otus.dt <- as.data.frame(some_ps@otu_table@.Data)
  jacco.mean <- as.list(vegdist(otus.dt, method = distance))
  return(jacco.mean)
}


bray.list <- sapply(seq.list, bray.steppo.nm)
bray.list <- bray.steppo.nm(0, ps = some_ps)
bray.list
length(bray.list)

wuni.list <- sapply(seq.list, wuni.steppo)
uni.list <- sapply(seq.list, uni.steppo)
bray.list
plot(bray.list)
plot(wuni.list)
plot(uni.list)

vegdist(otus.dt, method = "bray") %>% as.list

# train taxa reference set

getwd()
library(tidyverse)
library(DECIPHER)
seqs_path <- "~/storage/somebases/SILVA_138_SSURef_tax_silva.fasta"

seqs <- readDNAStringSet(seqs_path)
head(seqs)
taxid <- NULL
groups <- names(seqs)
head(groups.2)
groups.2 <- gsub("(.*\\..*\\.[0-9]+ )(.*)", "\\2", groups)
groups.2
groups[922]
groupCounts <- table(groups.2)
u_groups <- names(groupCounts)
length(u_groups)
length(groups.2)

# try some oldfasion methods from the 1000th tutorial
library(phyloseq)
library(DESeq2)

ps.s.bf.57@sam_data




