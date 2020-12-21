''' use script CRT_Functions_v1.1.R from https://github.com/ShadeLab/ConditionallyRareTaxa.git

see also https://github.com/rachelkaminsky2691/RK_CRT.git
'''

library(ggplot2)
library(ggthemes)
library(vegan)
library(phyloseq)
library(dplyr)

library(TSA)
source("CRT_Functions_v1.1.R")
SimpleRareToPrev.f(otu_fp="av_tp_r_otu.txt", abund_thresh=0.01, abund_thresh_ALL=FALSE, b_thresh=0.90, rdp_lastcol=FALSE)


write.table(av_tp_r_dat, "av_tp_r_otu.txt", quote = FALSE, sep = "\t")

av_tp_r_dat <- data.frame(t(otu_table(AY.ps)))
write.table(av_tp_r_dat, "AY.ps.txt", quote = FALSE, sep = "\t")
SimpleRareToPrev.f(otu_fp="AY.ps.txt", abund_thresh=0.001, abund_thresh_ALL=FALSE, b_thresh=0.90, rdp_lastcol=FALSE)

'''
Copy OTU column from results file, paste into a new spreadsheet and save as a csv

Prune original OTU table to only include CRT, save to directory 

Create csv file with all variable crt columns
'''
library(data.table)

crt.all <-as.data.frame(fread("CRT.csv"))
C.ps.with.crt <- prune_taxa(crt.all$C, AC.ps) # Exactly only CRT contain

pop_taxa = function(physeq, badTaxa){
    allTaxa = taxa_names(physeq)
    myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
    return(prune_taxa(myTaxa, physeq))
}

library(ape)
library(vegan)

AY.ps.without.crt <- pop_taxa(AY.ps, crt.all$AY)

dist.AY.ps.with.crt <- distance(AY.ps.with.crt, "bray")
dist.AY.ps.without.crt <- distance(AY.ps.without.crt, "bray")

mantel(dist.AY.ps.with.crt, dist.AY.ps.without.crt, method="pearson")

AY.ps.with.crt.pcoa <- pcoa(dist.AY.ps.with.crt, correction="none")
AY.ps.without.crt.pcoa <- pcoa(dist.AY.ps.without.crt, correction="none")
p1 <- plot_ordination(AY.ps, AY.ps.without.crt.pcoa, type = "samples", axes = 1:2,color="Site")
p2 <- plot_ordination(AY.ps, AY.ps.with.crt.pcoa, type = "samples", axes = 1:2,color="Site")
mantel(dist.AY.ps.with.crt, dist.AY.ps.without.crt, method="pearson")

library(ggpubr)

p <- ggarrange(p2, p1,  labels= c("CRT","все остальное"), ncol = 2, nrow = 1,label.y = 0.08)