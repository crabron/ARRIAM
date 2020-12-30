#!/usr/bin/env Rscript

library("Deseq2")





Des.Phylo <-  function(ps, Taxa){
    ps <- taxa_level(ps, Taxa)
    diagdds <- phyloseq_to_deseq2(ps, ~ Repeats)
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")    
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
   