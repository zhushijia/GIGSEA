setwd('/Users/shijiazhu/Dropbox/mypaper/GIGSEA/codes/data')

pathway <- readLines('KEGG/c2.cp.kegg.v4.0.symbols.gmt')
tf <- readLines('TF/c3.tft.v4.0.symbols.gmt')
mirna <- readLines('miRNA/c3.mir.v4.0.symbols.gmt')

MSigDB.KEGG.Pathway.gs <- buildGeneset( pathway )
MSigDB.TF.gs <- buildGeneset( tf )
MSigDB.miRNA.gs <- buildGeneset( mirna )


MSigDB.KEGG.Pathway.net = geneSet2sparseMatrix( term=MSigDB.KEGG.Pathway.gs$term , geneset=MSigDB.KEGG.Pathway.gs$geneset )
MSigDB.TF.net = geneSet2sparseMatrix( term=MSigDB.TF.gs$term , geneset=MSigDB.TF.gs$geneset )
MSigDB.miRNA.net = geneSet2sparseMatrix( term=MSigDB.miRNA.gs$term , geneset=MSigDB.miRNA.gs$geneset )

MSigDB.KEGG.Pathway.annot = MSigDB.KEGG.Pathway.gs[,1:2]
MSigDB.KEGG.Pathway.annot$totalGenes = as.integer(colSums(MSigDB.KEGG.Pathway.net))
MSigDB.TF.annot = MSigDB.TF.gs[,1:2]
MSigDB.TF.annot$totalGenes = as.integer(colSums(MSigDB.TF.net))
MSigDB.miRNA.annot = MSigDB.miRNA.gs[,1:2]
MSigDB.miRNA.annot$totalGenes = as.integer(colSums(MSigDB.miRNA.net))

MSigDB.KEGG.Pathway = list(net=MSigDB.KEGG.Pathway.net,annot=MSigDB.KEGG.Pathway.annot)
MSigDB.TF = list(net=MSigDB.TF.net,annot=MSigDB.TF.annot)
MSigDB.miRNA = list(net=MSigDB.miRNA.net,annot=MSigDB.miRNA.annot)

setwd('/Users/shijiazhu/Dropbox/mypaper/GIGSEA/codes/GIGSEA')
save(MSigDB.KEGG.Pathway,file='data/MSigDB.KEGG.Pathway.rda')
save(MSigDB.TF,file='data/MSigDB.TF.rda')
save(MSigDB.miRNA,file='data/MSigDB.miRNA.rda')

