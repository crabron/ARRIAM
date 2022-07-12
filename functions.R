
# установить ampvis2(ставится с гитхаба, в CRAN и Bioconductor вроде нет), tidyverse -- https://www.tidyverse.org/ там много классных штук
# https://stepik.org/course/724/promo - про dplyr - уроки 1.5 - 1.6 оч было бы неплохо если бы ты их прошел(без фанатизма)
# https://www.rstudio.com/resources/cheatsheets/ - и вот это - распечатываешь(в цвете!) про те пакеты которыми ты пользуешься и тогда становится понятно про что гуглить

phyloseq_to_amp <- function(ps){
  require(ampvis2)
  require(tibble)
  require(phyloseq)
  colnames(ps@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
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

# прогони то же, но с разряжением(разницы быть не должно, но ты попробуй)
plot_alpha_w_toc <- function(ps, group, metric) {
  
  require(phyloseq)
  require(ggplot2)
  
  ps_a <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  er <- estimate_richness(ps_a)
  df_er <- cbind(ps_a@sam_data, er)
  df_er <- df_er %>% select(c(group, metric))
  stat.test <- aov(as.formula(paste0(metric, "~", group)), data = df_er) %>%
    rstatix::tukey_hsd()
  y <- seq(max(er[[metric]]), length=length(stat.test$p.adj.signif[stat.test$p.adj.signif != "ns"]), by=max(er[[metric]]/20))
  
  plot_richness(ps_a, x=group, measures=metric) + 
    geom_boxplot() +
    geom_point(size=1.2, alpha=0.3) +
    stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", 
      y.position = y,
      hide.ns=TRUE) +
    theme_light() + 
    scale_color_brewer(palette="Dark2") +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank()) +
    labs(y=paste(metric, "index")) 
}

beta_custom_norm_NMDS_elli_w <- function(ps, seed = 7888, normtype="vst", Color="What", Group="Repeat"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  library(ggforce)
  
  
  # попробуй заменить NMDS на PCoA
  # https://joey711.github.io/phyloseq/plot_ordination-examples.html и вот здесь куча вариантов ординаций в phyloseq - попробуй нарисовать отдельно для таксонов, а не образцов(не моей функцией, у меня зависнет, и лучше PCoA - так быстрее))
  ordination.b <- ordinate(ps, "NMDS", "bray")
  mds <- as.data.frame(ordination.b$points)
  p  <-  plot_ordination(ps,
                         ordination.b,
                         type="sample",
                         color = Color,
                         # title="NMDS - Bray-Curtis",
                         title=NULL,
                         axes = c(1,2) ) + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    geom_point(size = 3) +
    annotate("text",
             x=min(mds$MDS1) + abs(min(mds$MDS1))/10,
             y=max(mds$MDS2),
             label=paste0("Stress -- ", round(ordination.b$stress, 3))) +
    geom_mark_ellipse(aes_string(group = Group, label = Group),
                      label.fontsize = 10,
                      label.buffer = unit(2, "mm"),
                      label.minwidth = unit(5, "mm"),
                      con.cap = unit(0.1, "mm"),
                      con.colour='gray') +
    theme(legend.position = "none") 
  
  return(p)
}

# бета с разряжением -- нужен seed! 

beta_custom_norm_NMDS_elli_rar <- function(ps, seed = 7888, normtype="vst", Color="What", Group="Repeat"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  library(ggforce)
  
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  #variant with only bray
  
  ps.rar <- rarefy_even_depth(ps)
  
  ordination.b <- ordinate(ps.rar, "NMDS", "bray")
  p  <-  plot_ordination(ps.rar,
                         ordination.b,
                         type="sample",
                         color = Color,
                         # title="NMDS - Bray-Curtis",
                         title=NULL,
                         axes = c(1,2) ) + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    geom_point(size = 3) +
    geom_mark_ellipse(aes_string(group = Group, label = Group),
                      label.fontsize = 10,
                      label.buffer = unit(2, "mm"),
                      label.minwidth = unit(5, "mm"),
                      con.cap = unit(0.1, "mm"),
                      con.colour='gray') +
    theme(legend.position = "none") 
  
  return(p)
}


# Variance stabilisation из пакета DESeq2 

beta_custom_norm_NMDS_elli <- function(ps, seed = 7888, normtype="vst", Color="What", Group="Repeat"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  library(ggforce)
  
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab
  #variant with only bray
  
  if (exists("ps.varstab") == FALSE){
    diagdds = phyloseq_to_deseq2(ps, ~ Site)                  
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
  }
  
  ordination.b <- ordinate(ps.varstab, "NMDS", "bray")
  p  <-  plot_ordination(ps.varstab,
                         ordination.b,
                         type="sample",
                         color = Color,
                         # title="NMDS - Bray-Curtis",
                         title=NULL,
                         axes = c(1,2) ) + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    geom_point(size = 3) +
    geom_mark_ellipse(aes_string(group = Group, label = Group),
                      label.fontsize = 10,
                      label.buffer = unit(2, "mm"),
                      label.minwidth = unit(5, "mm"),
                      con.cap = unit(0.1, "mm"),
                      con.colour='gray') +
    theme(legend.position = "none") 
  
  return(p)
}


plot_rich_reads_samlenames_lm <- function(physeq){
  rish <- estimate_richness(physeq, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(physeq))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Repeat"] <-unlist(purrr::map(stringr::str_split(rownames(physeq@sam_data), "straw-its2-", 2), function(x) x[[2]]))
  reads.summary["Site"] <- physeq@sam_data$Day
  library(ggrepel)
  require(ggforce)
  p1 <- ggplot(data=reads.summary) + 
    geom_point(aes(y=otus, x=log2(reads), color=Site),size=3) + 
    geom_text_repel(aes(y=otus, x=log2(reads), label=paste0(Site, "_", Repeat))) + 
    theme_bw() +
    geom_smooth(aes(y=otus, x=log2(reads), fill=Site, color=Site),method=lm, se=FALSE, ymin = 1) + 
    scale_x_continuous(sec.axis = sec_axis(sec.axis ~ 2**.)) 

  return(p1)
}

# рабочая функция, рисует барплоты, можешь посмотреть как филы отличаются от образца к образцу
bargraph <- function(ps, rank, threshold=0.05){
  require(dplyr)
  require(ggplot2)
  require(phyloseq)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps2 <- tax_glom(ps, taxrank = rank)
  ps3 = transform_sample_counts(ps2, function(x) x / sum(x) )
  data <- psmelt(ps3) # create dataframe from phyloseq object
  data$Plot <- as.character(data[,rank]) # convert to character
  data$Plot[data$Abundance < threshold] <- paste0("<", threshold, " abund.")
  medians <- data %>% group_by(Plot) %>% mutate(median=median(data$Abundance))
  remainder <- medians[medians$median <= threshold,]$Plot
  
  # create palette long enough for our data
  base.palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", 
                    "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", 
                    "darksalmon", "dodgerblue3", "steelblue1", "darkgoldenrod1", "brown1", "cyan1", "darkgrey")
  required.colors <- nlevels(factor(data$Plot))
  repeats = required.colors %/% length(base.palette) + 1
  palette <- rep(base.palette, length.out = repeats * length(base.palette))
  
  p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Plot))
  p + geom_bar(aes(), stat="identity", position="stack") + theme_light() +
    scale_fill_manual(values = palette) +
    theme(legend.position="bottom") + guides() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle('ITS2 - Phylum level')
  
}


amp_heatmap(amp.cc, tax_show = 30, group_by = "Association", tax_aggregate = "Genus", tax_add = "Phylum", tax_class = "Proteobacteria") +
  theme_bw() + 
  theme(text = element_text(size=15), legend.position = "none") + 
  theme(axis.text.x=element_text(angle=45,hjust=1))




ps_object <- ps 
ps_object <- subset_taxa(ps_object, Phylum != "NA")

ps_object@tax_table[is.na(ps_object@tax_table)] <- TRUE
ps_object <- subset_taxa(ps_object,
                         !(Family  == "Mitochondria" |
                             Class   == "Chloroplast" |
                             Order   == "Chloroplast"))
ps_object@tax_table <- dplyr::na_if(ps_object@tax_table, TRUE)
ps.f <- ps_object



ps_f_ph_1 <- tax_glom(ps.f, taxrank="Genus")

plot_tree(ps_f_ph_1,
          ladderize="left",
          color="Phylum", 
          label.tips="Class",
          base.spacing=0.015,
          text.size = 5)

