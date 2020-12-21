ggtree(kim.tree.16) + geom_tiplab(align = TRUE) + xlim(0, 0.008) + geom_hilight(node=23, fill="darkgreen", alpha=0.4) +  geom_text(aes(label=node), hjust=-5, size = 2.6) + geom_balance(node=99, fill='steelblue', color='white', alpha=0.6, extend=1)


ggtree(kim.tree.16, branch.length = "none")+ geom_tiplab(hjust = -0.5) + xlim(0, 70) + geom_hilight(node=13, fill="darkgreen", alpha=0.2,  extend=-2) +  geom_text(aes(label=node), hjust=-5, size = 2.6)

 ggtree(kim.tree.16) + geom_tiplab(align = TRUE, hjust=-0.3) + xlim(0, 0.02) +  geom_text(aes(label=node), hjust=-5, size = 2.6) + 
+     geom_cladelabel(node=2, label="vav-dag", color="red2", align=TRUE) + geom_cladelabel(node=1, label="rlt", color="red2", align=TRUE) + geom_cladelabel(node=2, label="vav-dag", color="red2", align=TRUE) + geom_cladelabel(node=3:4, label="vav-no", color="red2", align=TRUE) + geom_cladelabel(node=5:21, label="rst", color="red2", align=TRUE)

View(fortify(kim.tree.16))

http://bioconductor.org/packages/release/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html#annotate-clades
