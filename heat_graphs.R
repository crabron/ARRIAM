library(qgraph)

col.d <- otus
colnames(col.d) <- NULL
cor.col.d <- as.data.frame(cor(col.d))
Q <- qgraph(cor.col.d, layout = "spring")
col.d <- cbind(colnames(cor.col.d), taxa)
write.table(col.d, file = "storage/mal/heatmepe/col.result.tsv", sep= "\t", col.names = NA, quote=FALSE)

otus.filt <- otus[, !(c("CP009247.1399865.1401405","CRKN01000054.4139.5702", "LQOW01000023.119.1666"))]

otus.filt <-  otus[, -which(otus %in% c("CP009247.1399865.1401405","CRKN01000054.4139.5702", "LQOW01000023.119.1666"))]

data <- data[, !(names(data) %in% delete)] 
otus.data <- as.data.frame(otus)
otus.filt <-   select(otus.data, -c("CP009247.1399865.1401405","CRKN01000054.4139.5702", "LQOW01000023.119.1666"))
colnames(otus.filt) <- NULL
Q <- qgraph(cor.col.d, layout = "spring")
col.d.filt <-  col.d[-which(otus %in% c("CP009247.1399865.1401405","CRKN01000054.4139.5702", "LQOW01000023.119.1666")),]
write.table(col.d.filt, file = "~/storage/mal/heatmepe/col.result.tsv", sep= "\t", col.names = NA, quote=FALSE)

Q <- qgraph(cor(otus.filt), layout = "spring")
col.d <- cbind(colnames(otus.filt), taxa)
col.d.filt <-  select(col.d, -c("CP009247.1399865.1401405","CRKN01000054.4139.5702", "LQOW01000023.119.1666"))
write.table(col.d.filt, file = "~/storage/mal/heatmepe/col.result.tsv", sep= "\t", col.names = NA, quote=FALSE)

taxa.new <- taxa[!rownames(taxa) %in% "CP009247.1399865.1401405", ]
taxa.new <- taxa.new[!rownames(taxa.new) %in% "CRKN01000054.4139.5702", ]
taxa.new <- taxa.new[!rownames(taxa.new) %in% "LQOW01000023.119.1666", ]
col.d <- cbind(colnames(otus.filt), taxa.new)

# another  work
#write fasta from dataframe
library("seqRFLP")
dataframe2fas(briefToSeq.df, file="rep.fasta")

write.table(ps@tax_table), file = "otu_table.tsv", sep= "\t", col.names = NA, quote=FALSE)
