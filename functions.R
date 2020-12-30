
veganifyOTU <- function(physeq){
  require(phyloseq)
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}


pop.taxa <- function(physeq, badTaxa){
  require(phyloseq)
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}


delete_mit_chl <- function(ps){
  require(phyloseq)
  badTaxa <- taxa_names(subset_taxa(ps, Order=="Chloroplast"))
  ps <- pop.taxa(ps, badTaxa)
  badTaxa <- taxa_names(subset_taxa(ps, Family=="Mitochondria"))
  ps <- pop.taxa(ps, badTaxa)
  return(ps)
}

delete_mit_chl_for_rdp <- function(ps){
  require(phyloseq)
  badTaxa <- taxa_names(subset_taxa(ps, Class=="Chloroplast"))
  ps <- pop.taxa(ps, badTaxa)
  return(ps)
}

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


kate.ggcca.species <- function(vare.cca){
  
  require(ggvegan)
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.4,
                                                                                                                                                                                                         color = "blue",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
  return(p)
}


phyloseq_to_amp <- function(ps){
    require(ampvis2)
    require(tibble)
    require(phyloseq)
    OTU1 = as(otu_table(ps), "matrix")
    OTU1 <- t(OTU1)
    OTUdf = as.data.frame(OTU1)
    taxa.ps <- as(tax_table(ps), "matrix")
    taxa.df = as.data.frame(taxa.ps)
    my_otu_table <- merge(OTUdf, taxa.df, by=0)
    my_otu_table <- column_to_rownames(my_otu_table, var="Row.names")

    my_metadata <- as_tibble(sample_data(ps), rownames=NA)
    my_metadata <- rownames_to_column(my_metadata,var = "SampleID")
    my_tree <- phy_tree(ps)
    amp.ps <- amp_load(otutable = my_otu_table, metadata = my_metadata, tree = my_tree)
    return(amp.ps)
}

phyloseq_to_amp_without_tree <- function(ps){
  require(ampvis2)
  require(tibble)
  require(phyloseq)
  OTU1 = as(otu_table(ps), "matrix")
  OTU1 <- t(OTU1)
  OTUdf = as.data.frame(OTU1)
  taxa.ps <- as(tax_table(ps), "matrix")
  taxa.df = as.data.frame(taxa.ps)
  my_otu_table <- merge(OTUdf, taxa.df, by=0)
  my_otu_table <- column_to_rownames(my_otu_table, var="Row.names")
  
  my_metadata <- as_tibble(sample_data(ps), rownames=NA)
  my_metadata <- rownames_to_column(my_metadata,var = "SampleID")
  amp.ps <- amp_load(otutable = my_otu_table, metadata = my_metadata)
  return(amp.ps)
}


beta_for_Al_with_norm_NMDS <- function(ps, seed = 5433){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)

  
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Description)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
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
  # first row for 1-2 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Al", shape="Inoculation", title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Al", shape="Inoculation", title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Al", shape="Inoculation", title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  

  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1)
  
  return(p.all)
}


beta_for_Al_with_NMDS <- function(ps5){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)

  
  #beta and ordination
  
  ordination.b <- ordinate(ps, "NMDS", "bray")
  ordination.u <- ordinate(ps, "NMDS", "unifrac")
  ordination.w <- ordinate(ps, "NMDS", "wunifrac")
  
  #plotting
  # first row for 1-2 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Al", shape="Inoculation", title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Al", shape="Inoculation", title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Al", shape="Inoculation", title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1)
  
  return(p.all)
}

permanova.inoculation <- function(ps, dist = "bray"){
  require(phyloseq)
  require(vegan)
  dist <- distance(ps, dist)
  metadata <- as(sample_data(ps), "data.frame")
  ad <- adonis2(dist ~ Inoculation, data = metadata)
  return(ad)
}


permanova.custom <- function(ps, split = FALSE, factor = "Type", factor_variance = "pos", dist = "bray", formula = dist ~ Repeats){
  require(phyloseq)
  require(vegan)
  while (split){
    ps.Cd <- prune_samples(sample_data(ps)$Cd %in% c("pos"), ps)
    ps.Cd <- prune_taxa(taxa_sums(ps.Cd) > 0, ps.Cd)    
  }
  dist <- distance(ps, dist)
  metadata <- as(sample_data(ps), "data.frame")
  ad <- adonis2(dist ~ split_factor, data = metadata)
  return(ad)
}


permanova.al <- function(ps, dist = "bray"){
  require(phyloseq)
  require(vegan)
  dist <- distance(ps, dist)
  metadata <- as(sample_data(ps), "data.frame")
  ad <- adonis2(dist ~ Al, data = metadata)
  return(ad)
}


al.ggcca.sites <- function(vare.cca){
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

beta_for_Al_with_alter_norm_NMDS <- function(ps, seed = 235){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  
  
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds <- phyloseq_to_deseq2(ps, ~ Description)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
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
  # first row for 1-2 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Al", shape="Inoculation", title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Al", shape="Inoculation", title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Al", shape="Inoculation", title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1)
  
  return(p.all)
}


al.ggcca.species <- function(vare.cca){
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.8,
                                                                                                                                                                                                             color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
}

Des.Al <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ Drought)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Drought), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
}      

Des.Tax = function(ps, Taxa){
  ps <- taxa_level(ps, Taxa)
  diagdds = phyloseq_to_deseq2(ps, ~ Drought)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Drought), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df, as.data.frame(aggdata)[rownames(res.df),])
  return(nice)
}  

mantel.wrap.veg <- function(ps){}

beta_for_Plates_with_norm_NMDS <- function(ps, seed = 5433){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  
  
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
#  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  
  ps.rand <- rarefy_even_depth(ps, rngseed = seed)
  
  #beta and ordination
  
  ordination.b <- ordinate(ps.varstab, "NMDS", "bray")
  ordination.u <- ordinate(ps.rand, "NMDS", "unifrac")
  ordination.w <- ordinate(ps.varstab, "NMDS", "wunifrac")
  
  #plotting
  # first row for 1-2 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Type", shape="Soil", title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Type", shape="Soil", title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Type", shape="Soil", title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1)
  
  return(p.all)
}

beta_for_Plates <- function(ps){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  
  #beta and ordination
  
  ordination.b <- ordinate(ps, "PCoA", "bray")
  ordination.u <- ordinate(ps, "PCoA", "unifrac")
  ordination.w <- ordinate(ps, "PCoA", "wunifrac")
  
  #plotting
  # first row for 1-2 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Type", shape="Soil", title="PCoA - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Type", shape="Soil", title="PCoA - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Type", shape="Soil", title="PCoA - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #second row for 1-3 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Type", shape="Soil", title="PCoA - Bray", 
                      axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p2.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample",color="Type", shape="Soil", title="PCoA - unifrac", 
                      axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p2.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Type", shape="Soil", title="PCoA - wunifrac", 
                      axes = c(1,3) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p2.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, p2.1, p2.2, p2.3, ncol = 3 , nrow = 2)
  
  return(p.all)
}

beta_for_Svir_with_norm_NMDS <- function(ps, seed = 5433){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  
  
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  #  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  
  ps.rand <- rarefy_even_depth(ps, rngseed = seed)
  
  #beta and ordination
  
  ordination.b <- ordinate(ps.varstab, "NMDS", "bray")
  ordination.u <- ordinate(ps.rand, "NMDS", "unifrac")
  ordination.w <- ordinate(ps.varstab, "NMDS", "wunifrac")
  
  #plotting
  # first row for 1-2 axis
  
  p = plot_ordination(ps, ordination.b, type="sample", color="Type", shape="Soil", title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color="Type", shape="Soil", title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color="Type", shape="Soil", title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1)
  
  return(p.all)
}

beta_for_Plates_with_norm_NMDS <- function(ps, seed = 5433){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  
  
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  #  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  
  ps.rand <- rarefy_even_depth(ps, rngseed = seed)
  
  #beta and ordination
  
  ordination.b <- ordinate(ps.varstab, "NMDS", "bray")
  ordination.u <- ordinate(ps.rand, "NMDS", "unifrac")
  ordination.w <- ordinate(ps.varstab, "NMDS", "wunifrac")
  
  #plotting
  # first row for 1-2 axis
  
  p = plot_ordination(ps.varstab, ordination.b, type="sample", color="Type", shape="Soil", title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps.rand, ordination.u, type="sample", color="Type", shape="Soil", title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps.varstab, ordination.w, type="sample", color="Type", shape="Soil", title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1)
  
  return(p.all)
}

beta_custom_norm_NMDS <- function(ps, seed = 6788, normtype="vst", color="Cd", shape="Drought"){
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
  p = plot_ordination(ps, ordination.b, type="sample", color, shape, title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color, shape, title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color, shape, title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1, common.legend = TRUE, legend = "right")
  
  return(p.all)
}
