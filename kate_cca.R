#офигительный cca-biplot со стрелочками в ggplot2

library(vegan)
library(ggplot)
library(dplyr)
library(ggrepel)

fdat <- fortify(vare.cca)
p.species <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + geom_point(mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + 
  geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=3) 
p.species + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.2, color = "red",arrow = arrow(angle = 3))  + 
  theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))




little.cca <- cca(otus.ps.vegan ~ pHh2o , data=metadata)
mod <- ordistep(little.cca, scope = formula(vare.norm.cca))




#как построить cca в vegan из phyloseq

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
map <- as.data.frame(t(map))

vare.cca <- cca(otus.ps.vegan ~ Depth + Lake_dist + N_NO3 + HA + EA + pHh2o + Carbon , data=map)


#рисовалка биплота базовыми средствами
plot(vare.cca,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,4),scaling=2)
points(vare.cca,disp="sites",pch=21,col="red",bg="red",cex=1)
text(vare.cca,"sites",pos=4,axis.bp=TRUE)

#ссылочки

https://www.researchgate.net/post/Can_any_one_help_me_with_the_interpretation_of_CCA_plot
http://cc.oulu.fi/~jarioksa/opetus/metodi/mmmbeam2.pdf


fdat <- fortify(vare.cca)
p.sites <- ggplot(fdat %>% filter(Score %in% c("sites","biplot"))) + geom_point(mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=3) 
p.sites + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.2, color = "red",arrow = arrow(angle = 3))  + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50")))



vare.cca
anova(vare.cca, by="terms")
anova(vare.cca, by="mar")
little.cca <- cca(otu ~pHh2o , data=mapping)
mod <- ordistep(little.cca, scope = formula(vare.cca))
mod$anova
vif.ccaa$vare.cca

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


p.sites <- ggplot(fdat_ch %>% filter(Score %in% c("sites","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "sites"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat_ch %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=1, color = "black",arrow = arrow(angle = 3))  + 
  geom_label(data = fdat_ch %>% dplyr::filter(Score == "biplot"), aes(x= CCA1, y=CCA2 + 0.0, label = tax)) +  geom_text_repel(data=fdat_ch %>% filter(Score == "sites"),mapping=aes(x=CCA1, y=CCA2, label= Label),size=4)
p <-  p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
p

kate.ggcca.species <- function(vare.cca){
  
  require(ggvegan)
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + 
    geom_point(data = fdat %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=1,  geom_label(data = fdat_ch %>% dplyr::filter(Score == "biplot"), aes(x= CCA1, y=CCA2, label = tax)) + 
    color = "blue",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
  return(p)
}


load("~/storage/r-base_files/svir.RData")

fdat <- fortify(vare.cca)
kate.ggcca.species(vare.cca.wd)
write.table(fdat, file = "fdat.txt", sep= "\t", col.names = NA, quote=FALSE)

mapp <- read.csv("map_with_dist_depth.csv" , header=TRUE, sep="\t")
map <- data.frame(row.names="ID", mapp)
View(map)

fdat_ch <- read.csv("fdat.txt" , header=TRUE, sep="\t")

p.sites <- ggplot(fdat_ch %>% filter(Score %in% c("species","biplot"))) + geom_point(data = fdat_ch %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2)) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=1,
                                                                                                                                                                                                           color = "black",arrow = arrow(angle = 3))  + geom_label(data = fdat_ch %>% dplyr::filter(Score == "biplot"), aes(x= CCA1, y=CCA2, label = tax)) +
  geom_text_repel(data = fdat_ch %>% dplyr::filter(Score == "species"), aes(x=CCA1, y=CCA2, label= tax, colour=phylum),size=4) 
p <- p.sites + theme(legend.position = "right", legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"))
p
library(ggplot2)
p2 <- p + scale_color_manual(values=c("#5EA09E", "#6B440B", "#F499C2", "#7DA9D8", "#FEF898", "#0000FE", "#F27204", 
                                "#FE0000", "#FEFF00", "#02F30E", "#F89679", "#7F8000", "black","black"))
p
load("~/storage/r-base_files/mal_cca.RData")
saveRDS(p2, file = "~/storage/see_some_results/kate_cca_plot")
                     
adonis.pair<-function(dist.mat,Factor,nper=1000,corr.method="fdr"){
  require(vegan)
  as.factor(Factor)
  comb.fact<-combn(levels(Factor),2)
  pv<-NULL
  R2<-NULL
  SS<-NULL
  MeanSqs<-NULL
  F.Model<-NULL
  for (i in 1:dim(comb.fact)[2]){
    model.temp<-adonis(as.dist(as.matrix(dist.mat)[Factor==comb.fact[1,i] | Factor==comb.fact[2,i],Factor==comb.fact[1,i] | Factor==comb.fact[2,i]])~Factor[Factor==comb.fact[1,i] | Factor==comb.fact[2,i]],permutations=nper)
    pv<-c(pv,model.temp$aov.tab[[6]][1])
    R2<-c(R2,model.temp$aov.tab$R2[1])
    SS<-c(SS,model.temp$aov.tab[[2]][1])
    MeanSqs<-c(MeanSqs,model.temp$aov.tab[[3]][1])
    F.Model<-c(F.Model,model.temp$aov.tab[[4]][1])}
  pv.corr<-p.adjust(pv,method=corr.method)
  dt <- data.frame(combination=paste(comb.fact[1,],comb.fact[2,],sep=" <-> "),SumsOfSqs=SS,MeanSqs=MeanSqs,F.Model=F.Model,R2=R2,P.value=pv,P.value.corrected=pv.corr)}
  return(dt)
}

ps.dist <- vegdist(otus.ps.vegan, "bray")
Factor <- ps@sam_data$Repeats
adonis.pair(ps.dist, ps@sam_data$Repeats)
comb.fact<-combn(levels(Factor),2)
