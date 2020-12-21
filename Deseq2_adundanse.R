diagdds_more = function(ps){
    diagdds = phyloseq_to_deseq2(ps, ~ Repeats)
    diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
    res = results(diagdds)
    res = res[order(res$padj, na.last=NA), ]
    sigtab = res[(res$padj < 0.1), ]
    sigtab = sigtab[(sigtab$log2FoldChange > 1),]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
    return(sigtab)
}


diagdds_less = function(ps){
    diagdds = phyloseq_to_deseq2(ps, ~ Repeats)
    diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
    res = results(diagdds)
    res = res[order(res$padj, na.last=NA), ]
    sigtab = res[(res$padj < 0.1), ]
    sigtab <- sigtab[(sigtab$log2FoldChange < 1),]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
    return(sigtab)
}

change.prop <-prop.table(table(sigtab$Phylum))

boxplot(log10(assays(diagdds)[["cooks"]]), range=0, las=2)

cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
    
    
    
cook <- dds@assays[["cooks"]]
                  
For Lise:
                  

diagdds_Lise = function(ps, name){
    diagdds <- phyloseq_to_deseq2(ps, ~ Repeats)
    
    samp <-sample_data(ps)
    dds.counts <- diagdds@assays@.xData$data$counts
    dds.counts.df <- as.data.frame(dds.counts)
    aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Repeats), median))
    colnames(aggdata) <- aggdata[1,]
    aggdata <- aggdata[-1,]
    res = results(diagdds)
    res.df <- as.data.frame(res)
    nice <- cbind(res.df, as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])
    
    return(nice)
}  

diagg.var <- diagdds_Lise(ps.var, "site1_diagdds.csv")                  
diagg.art <- diagdds_Lise(ps.art, "site2_diagdds.csv")                   
diagg.pse <- diagdds_Lise(ps.pse, "site3_diagdds.csv")                  
diagg.sph <-diagdds_Lise(ps.sph, "site4_diagdds.csv")                  
diagg.bac <-diagdds_Lise(ps.bac "site5_diagdds.csv")                   

> View(diagg.var)
> View(diagg.art)
> View(diagg.pse)
> View(diagg.sph)
> View(diagg.bac)                  
                  
nice <- cbind(as.data.frame(sigtab), as.data.frame(tax_table(ps.1)[rownames(sigtab),]), as.data.frame(aggdata[rownames(sigtab),])) 
                  
                  EF517956.1.1666
ggplot(data=depth.mut, aes(log(Conc_by_RealTime), ratio)) + geom_point() + ggrepel::geom_text_repel(data=subset(depth.mut, log(depth.mut$Conc_by_RealTime) < 16 ), aes(label=ID), size = 3)
                  
cooks.clean <- t(log10(assays(diagdds.ps.all.clean)[["cooks"]]))
cooks.clean <- rowMeans(cooks.clean, na.rm = TRUE)

diagdds_taxas = function(ps, taxa_level){
    physeq <- taxa_level(ps, taxa_level)
    diagdds <- phyloseq_to_deseq2(physeq, ~ Repeats)
    diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")
    res = results(diagdds)
    res.df <- as.data.frame(res)
    return(res.df)
}                    
                  
Des.Lise <- function(ps){
    otus.ps.vegan <- veganifyOTU(ps)
    metadata <- as(sample_data(ps), "data.frame")
    sim <- with(metadata, simper(otus.ps.vegan, Description))
    simper <- cbind(sim$s1903_In_s1903_In_Al$species,sim$s1903_In_s1903_In_Al$average)
    colnames(simper) <- c("ID","ave_sim")
    simper <- as.data.frame(simper, row.names = "ID")
    simper <- column_to_rownames(simper, var = "ID")
    diagdds <- phyloseq_to_deseq2(ps, ~ Description)
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")
    samp <-sample_data(ps)
    dds.counts <- diagdds@assays@.xData$data$counts
    dds.counts.df <- as.data.frame(dds.counts)
    aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Description), median))
    colnames(aggdata) <- aggdata[1,]
    aggdata <- aggdata[-1,]
    res = results(diagdds)
    res.df <- as.data.frame(res)
    nice <- cbind(res.df, simper[rownames(res.df),], as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])
    
    return(nice)
}                  
                  
 
Des.Norm <- function(ps){
    diagdds <- phyloseq_to_deseq2(ps, ~ Repeats)
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")
    dds.counts <- diagdds@assays@.xData$data$counts
    dds.counts.df <- as.matrix(dds.counts)
    ps.norm.dec <- phyloseq(otu_table(t(dds.counts.df), taxa_are_rows=FALSE), 
        sample_data(ps@sam_data), 
        tax_table(ps@tax_table@.Data),
        phy_tree(ps@phy_tree))
    ps.norm.dec
    nice <- ps.norm.dec
    
    return(nice)
}     
    
Des.Tax = function(ps, Taxa){
    ps <- taxa_level(ps, Taxa)
    diagdds <- phyloseq_to_deseq2(ps, ~ Description)
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")    
    samp <-sample_data(ps)
    dds.counts <- diagdds@assays@.xData$data$counts
    dds.counts.df <- as.data.frame(dds.counts)
    aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Description), median))
    colnames(aggdata) <- aggdata[1,]
    aggdata <- aggdata[-1,]
    res = results(diagdds)
    res.df <- as.data.frame(res)
    nice <- cbind(res.df, as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])
    
    return(nice)
}   
                  
Des.Phylo <-  function(ps, Taxa){
    ps <- taxa_level(ps, Taxa)
    diagdds <- phyloseq_to_deseq2(ps, ~ Description)
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")    
    samp <-sample_data(ps)
    dds.counts <- diagdds@assays@.xData$data$counts
    dds.counts.df <- as.data.frame(dds.counts)
    aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Description), median))
    colnames(aggdata) <- aggdata[1,]
    aggdata <- aggdata[-1,]
    res = results(diagdds)
    res.df <- as.data.frame(res)
    nice <- cbind(res.df, as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])
    
    return(nice)
}   

                  
                  
ps.1903 <- prune_taxa(taxa_sums(ps.1903) > 0, ps.1903)                  
diagddsraw = phyloseq_to_deseq2(ps.1903, ~ Description)                  
iagdds = estimateSizeFactors(diagddsraw, type="poscounts")                  
GPdds = estimateDispersions(iagdds, fitType = "local")                  
otu_table(ps.1903.varstab) <- otu_table(t(getVarianceStabilizedData(GPdds)), taxa_are_rows = FALSE)                  
  

ps.1903.mod <- prune_taxa(taxa_sums(ps.1903) > 10, ps.1903)              
diagdds <- phyloseq_to_deseq2(ps.1903.mod, ~ Description)                  
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
pst <- varianceStabilizingTransformation(diagdds, fitType="mean")               
pst.dimmed <- t(as.matrix(assay(pst))) 
pst.dimmed[pst.dimmed < 0.0] <- 0.0 
ps.varstab.mod <- ps.1903.mod            
otu_table(ps.varstab.mod) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)   
ordination.b <- ordinate(ps.varstab.mod, "PCoA", "bray")                  
p <- plot_ordination(ps.varstab.mod, ordination.b, type="sample", color="Al", shape="Inoculation", title="PCoA - Bray", axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
p + stat_ellipse( type="norm", alpha=0.7)              
                                    
                  
            
diagdds <- phyloseq_to_deseq2(ps.1903,  ~ Description)   
                  
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
pst <- varianceStabilizingTransformation(diagdds, fitType="mean")               
pst.dimmed <- t(as.matrix(assay(pst))) 
pst.dimmed[pst.dimmed < 0.0] <- 0.0 
ps.varstab.mod <- ps.1903.mod            
otu_table(ps.varstab.mod) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)                 
                  
Des.Al <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ Description)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
   diagdds = DESeq(diagdds)
    samp <-sample_data(ps)
    dds.counts <- diagdds@assays@.xData$data$counts
    dds.counts.df <- as.data.frame(dds.counts)
    aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Description), median))
    colnames(aggdata) <- aggdata[1,]
    aggdata <- aggdata[-1,]
    res = results(diagdds)
    res.df <- as.data.frame(res)
    nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
    return(nice)
}               
                  
  diagdds = phyloseq_to_deseq2(ps, ~ Description)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)              