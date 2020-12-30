setwd("~/storage/plates/plos/")
library(phyloseq)
ps <- readRDS(file = "~/storage/plates/ps.f.plos.rds") 
metadata <- ps@sam_data
metadata$Description <-  as.character(metadata$Repeats)
metadata$Description[metadata$Description == "getchinson_sod-podzolic"] <- "gSP"
metadata$Description[metadata$Description == "getchinson_black soil"] <- "gCZ"
metadata$Description[metadata$Description == "soil dna_sod-podzolic"] <- "bsSP"
metadata$Description[metadata$Description == "soil dna_black soil"] <- "bsCZ"
metadata$Description <- factor(metadata$Description, levels =  c("gSP", "gCZ", "bsSP", "bsCZ"))
sample_data(ps) <- metadata
library(DivNet)

ps.merged.repeats <- merge_samples(ps.f, "Description")
divnet_merged <-  divnet(ps.merged.repeats, ncores = 4)
divnet_merged$`bray-curtis`
divnet_phylum$shannon
divnet_asv <-  divnet(ps.f, ncores = 4)
divnet_asv
divnet_merged
saveRDS(divnet_merged, file = "~/storage/plates/divnet_merged.rds")
saveRDS(divnet_asv, file = "~/storage/plates/divnet_asv.rds")

p.little
beta.plot
p2
p.lm
