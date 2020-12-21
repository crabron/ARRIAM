library(phyloseq)
library(tibble)
library(dplyr)


col.mat <- rownames(as.data.frame(ps.all.w.clean@otu_table@.Data))
norm.prop.ps <- transform_sample_counts(ps.all.w.clean, function(OTU) OTU/sum(OTU))

s1 <- prune_samples(sample_data(ps.all.w.clean)$Site %in% c("siteI"), ps.all.w.clean)
s2 <- prune_samples(sample_data(ps.all.w.clean)$Site %in% c("siteII"), ps.all.w.clean)
s3 <- prune_samples(sample_data(ps.all.w.clean)$Site %in% c("siteIII"), ps.all.w.clean)
s4 <- prune_samples(sample_data(ps.all.w.clean)$Site %in% c("siteIV"), ps.all.w.clean)
s5 <- prune_samples(sample_data(ps.all.w.clean)$Site %in% c("siteV"), ps.all.w.clean)

s1 <- prune_taxa(taxa_sums(s1) > 5, s1)
s2 <- prune_taxa(taxa_sums(s2) > 5, s2)
s3 <- prune_taxa(taxa_sums(s3) > 5, s3)
s4 <- prune_taxa(taxa_sums(s4) > 5, s4)
s5 <- prune_taxa(taxa_sums(s5) > 5, s5)

mat.merged <- suppressMessages(full_join(as.data.frame(s1@otu_table@.Data),as.data.frame(s2@otu_table@.Data), name=TRUE))
mat.merged <- suppressMessages(full_join(mat.merged,as.data.frame(s3@otu_table@.Data), name=TRUE))
mat.merged <- suppressMessages(full_join(mat.merged,as.data.frame(s4@otu_table@.Data), name=TRUE))
mat.merged <- suppressMessages(full_join(mat.merged,as.data.frame(s5@otu_table@.Data), name=TRUE))
                                
rownames(mat.merged) <- col.mat
mat.merged.wna <- t(na.omit(t(mat.merged)))
col.filt <- colnames(mat.merged.wna)
filtered.ps <- prune_taxa(col.filt, norm.prop.ps)
amp.lise.out <- merge_samples(filtered.ps, "Repeats", fun=median)
lise.out <- cbind2(as.data.frame(t(amp.lise.out@otu_table@.Data)), as.data.frame(tax_table(ps.lise.out@tax_table@.Data)))
write.table(lise.out, file = "lise.out.csv", sep= "\t", col.names = NA, quote=FALSE)