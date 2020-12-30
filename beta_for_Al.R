beta_for_Al <- function(ps){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  
  #beta and ordination
  
  ordination.b <- ordinate(ps, "PCoA", "bray")
  ordination.u <- ordinate(ps, "PCoA", "unifrac")
  ordination.w <- ordinate(ps, "PCoA", "wunifrac")
  
  #plotting
  # first row for 1-2 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Al", shape="Inoculation", title="PCoA - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Al", shape="Inoculation", title="PCoA - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Al", shape="Inoculation", title="PCoA - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #second row for 1-3 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Al", shape="Inoculation", title="PCoA - Bray", 
                      axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p2.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Al", shape="Inoculation", title="PCoA - unifrac", 
                      axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p2.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Al", shape="Inoculation", title="PCoA - wunifrac", 
                      axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p2.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, p2.1, p2.2, p2.3, ncol = 3 , nrow = 2)
  
  return(p.all)
}