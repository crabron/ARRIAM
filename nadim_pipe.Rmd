library(phyloseq)
library(ggplot2)
library(ggpubr)

ordination <- ordinate(ps.f.norm,"PCoA", "unifrac")
p = plot_ordination(ps.f.norm, ordination, type="sample", color="Repeats", shape="Site", title="PCoA - unweighed Unifrac", axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 30)) + geom_point(size = 4) 
p1 <- p + stat_ellipse( type="norm", alpha=0.7)

ordination <- ordinate(ps.f.norm,"PCoA", "wunifrac")
p = plot_ordination(ps.f.norm, ordination, type="sample", color="Repeats", shape="Site", title="PCoA - Weighed Unifrac", axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 30)) + geom_point(size = 4) 
p2 <- p + stat_ellipse( type="norm", alpha=0.7)

ordination <- ordinate(ps.f.norm,"PCoA", "unifrac")
p = plot_ordination(ps.f.norm, ordination, type="sample", color="Repeats", shape="Site", title="PCoA - unweighed Unifrac", axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 30)) + geom_point(size = 4) 
p3 <- p + stat_ellipse( type="norm", alpha=0.7)

ordination <- ordinate(ps.f.norm,"PCoA", "wunifrac")
p = plot_ordination(ps.f.norm, ordination, type="sample", color="Repeats", shape="Site", title="PCoA - Weighed Unifrac", axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 30)) + geom_point(size = 4) 
p4 <- p + stat_ellipse( type="norm", alpha=0.7)

ggarrange(p1, p2, p3, p4, ncol = 2 , nrow = 2)


library(ggplot)
library(plyr)

lise.heat.1 <- read.csv("lise_heat_1.csv" , header=TRUE, sep="\t")

somedata <- ddply(lise.heat.1, .(Site),transform, pos = cumsum(Size) - (0.5 * Size))

ggplot() + geom_bar(aes(y = Size, x = Site, fill = OTU), data = somedata,  stat="identity", colour = "black") + geom_text(data = somedata, aes(x = Site, y = pos, label = X), size = 3) + scale_fill_gradient(low = "azure", high = "darkred", name = "Nitrososphaeraceae") +  theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y.right=element_blank(), axis.title.x=element_blank(),axis.text=element_text(size=13))