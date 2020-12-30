
#ps.pruned <- prune_taxa(taxa_sums(ps.clean) > 10, ps.clean ) 
#somedata.maj <- subset(somedata, rowSums(somedata) > 210)
#ps.merged <- merge_samples(ps.varstab, "Repeats", mean)
#otus.taxa.1[otus.taxa.1 == ""] <- NA
#lise.heat.1 <- read.csv("lise.heat.2.csv" , header=TRUE, sep="\t")
#lapply(otus.taxa$ID, lise.heat)

lise.heat <- function(x){
  require(ggplot2)
  require(plyr)
  taxafilename <- as.character(x)
  otu.vector <- otus.taxa[otus.taxa$ID == taxafilename, ]
  NonNAindex <- which(!is.na(otu.vector))
  lastNoNna <- max(NonNAindex)
  legend_name <- paste(as.character(otu.vector$Phylum), as.character(otu.vector[,lastNoNna]), sep = " -- ")
  somedata <- ddply(lise.heat.1, .(Site),transform, pos = cumsum(Size) - (0.5 * Size))
  #somedata$OTU <- rescale(somedata$OTU, to = c(0,1))
  otu.vector <- as.numeric(as.vector(otu.vector))[2:22]
  somedata$OTU <- otu.vector
  
  p <- ggplot() + geom_bar(aes(y = Size, x = Site, fill = OTU), data = somedata,  stat="identity",colour = "black", size = 0.3 ) + 
    geom_text(data = somedata, aes(x = Site, y = pos, label = Horizont), size = 3) + ggtitle(legend_name) +
    scale_fill_gradient2(low = "white",mid = "yellowgreen", high = "red4" , midpoint = max(somedata$OTU)/2) +
    theme(panel.border = element_blank(),panel.background = element_blank(),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y.right=element_blank(), 
    axis.title.x=element_blank(),axis.text=element_text(size=13)) 
    ggsave(taxafilename, p, device = "png", width = 20, height = 20, units = "cm")
  }