uni <- read.csv("~/storage/spades/PROKKA_11052019/uni.csv" , header=FALSE, sep="\t")
uni <- uni[-c(1),]
uni <- uni %>% remove_rownames %>% column_to_rownames(var="V1")


cog <- read.csv("~/storage/spades/PROKKA_11052019/COGforMINPATH.txt" , header=FALSE, sep="\t")





head(uni)
head(cog)


pv <- pathview(cpd.data = uni, pathway.id = cog$V2, species = "Streptomyces venezuelae", kegg.native = F)
View(cog)

pv <- pathview(cpd.data = uni, species = "Streptomyces venezuelae", kegg.native = F)

data(cpd.simtypes)
cpd.simtypes
data(gse16873.d)
View(gse16873.d)

ko <- read.csv("~/storage/spades/PROKKA_11052019/ko_table.csv" , header=FALSE, sep="\t")
ko <- ko[-c(1),]
ko <- ko %>% remove_rownames %>% column_to_rownames(var="V1")
View(ko)
data(demo.paths)
pv.out <- pathview(gene.data = ko , pathway.id = rownames(ko) ,species = "Streptomyces venezuelae", kegg.native = F)
sve