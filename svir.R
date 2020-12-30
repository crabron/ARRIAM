ps.f.i.pr.man <- prune_taxa(taxa_sums(ps.f.i) > 100, ps.f.i)

library(vegan)
some_ps <- ps.s.bf
ps.f.i.pr.man <- prune_taxa(taxa_sums(ps.f.i) > 100, ps.f.i)
otus.dt <- as.data.frame(some_ps@otu_table@.Data)
mean(vegdist(otus.dt, method = "jaccard"))

library(reticulate)
use_condaenv("pandas", conda = "/home/gladkov/.conda/envs/pandas")
use_python("/home/gladkov/.conda/envs/pandas/bin/python")

generator <- function(some_ps){
  min.reads <- min(taxa_sums(some_ps))
  new.min <- min.reads + 10
  some_ps <- prune_taxa(taxa_sums(some_ps) < new.min , some_ps) 
  if (min(taxa_sums(some_ps)) > 100)
    some_ps
  else
   NULL
}

ps.s.c <- prune_samples(sample_data(ps.s)$Horizont %in% c("C"), ps.s)
ps.s.c <- prune_taxa(taxa_sums(ps.s.c) > 0, ps.s.c) 

ps.s@sam_data$Horizont


seq.list <- seq(0, 500, by=10)
jac.steppo <- function(some_number){
  some_ps <- prune_taxa(taxa_sums(some_ps) > some_number , some_ps) 
  otus.dt <- as.data.frame(some_ps@otu_table@.Data)
  jacco.mean <- mean(vegdist(otus.dt, method = "jaccard"))
  return(jacco.mean)
}

some_ps <- ps.s
seq.list <- seq(0, 5000, by=10)
jacco.list <- sapply(seq.list, jac.steppo)
plot(jacco.list)

prune_taxa(taxa_sums(some_ps) > 5000 , some_ps)
