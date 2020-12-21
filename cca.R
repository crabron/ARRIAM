library(phyloseq)             
library(vegan)

#мержим по повторностям, ps.f - исходный phyloseq объект
ps.merged <- merge_samples(ps.f, "Description")
physeq <- ps.merged

# функция преобразует otu-table из phyloseq объекта в otu-table необходимый для пакета vegan
# естественно, если начальные данные не phyloseq объект, то этот кусок кода можно пропустить
veganifyOTU <- function(physeq){
  require(phyloseq)
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}


otus.ps.vegan <- veganifyOTU(physeq)
metadata <- as(sample_data(physeq), "data.frame")

#собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
#на вход - otu-table - matrix, metadata - data.frame
vare.cca <- vegan::cca(otus.ps.vegan ~  N + TOC + C + pH + K2O + P2O5, data=metadata)

#рисовалочка
kate.ggcca.sites <- function(vare.cca){
  require(ggvegan)
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  
  #вот этот участок кода рaботает на сторонней библиотеке ggvegan, довольно кривое решение, ниже вариант кода без использования ggvegan
  #
  fdat <- fortify(vare.cca)
  #
  
  fdat[fdat$Score == "sites",]
  
  p.sites <- ggplot(fdat %>% filter(Score %in% c("sites","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "sites"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + 
    geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.8, color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
}

#альтернативный кусок кода без использования ggrepel::fortify

require(tidyverse)

biplot <- as.data.frame(vare.cca$CCA$biplot)
wa <- as.data.frame(vare.cca$CCA$wa)

biplot <- rownames_to_column(biplot, "Label") %>% 
  add_column(Score = rep("biplot", length(rownames(biplot))))
wa <- rownames_to_column(wa, "Label") %>% 
  add_column(Score = rep("sites", length(rownames(wa))))
fdat_amazing <- rbind(biplot, wa)

# fdat_amazing - вместо fdat

kate.ggcca.sites(vare.cca)

# то же, но для организмов, а не сайтов

kate.ggcca.species <- function(vare.cca){
  
  require(ggvegan)
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + 
    geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.4, color = "blue",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
  return(p)
}

# для проверки модели использовать следующие функции:
# в целом 
anova(vare.cca)

#по отдельным осям
anova(vare.cca, by="terms")

#мультиколлениарность факторов(норма - значение менее 10)
vif.cca(vare.cca)

#для крайних значений модели
anova(vare.cca, by="mar")

#если слишком много факторов, то прогнать ordistep с разными алгоритмами