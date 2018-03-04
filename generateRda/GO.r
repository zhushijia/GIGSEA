library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)

go.term = toTable(GOTERM)
symbol2entrez = toTable(org.Hs.egSYMBOL)
go2entrez = toTable(org.Hs.egGO2EG)

uni_paste_GO = unique( paste( as.character(go.term[,1]),as.character(go.term[,3]) ,sep='_'  ) )
uni.go.term = t( sapply(strsplit(uni_paste_GO,'_'),function(x) { c(as.character(x[1]),as.character(x[2])) } ) )

range1 = match( as.character(go2entrez[,2]),as.character(uni.go.term[,1]) )
range2 = match( as.character(go2entrez[,1]),as.character(symbol2entrez[,1]) )
forGO = data.frame(symbol2entrez[range2,],go2entrez,uni.go.term[range1,])
annotGO = forGO[,c(2,4,8)]
colnames(annotGO) = c('symbol','go_id','go_term')

bp = toTable(GOBPOFFSPRING)
cc = toTable(GOCCOFFSPRING)
mf = toTable(GOMFOFFSPRING)
GOOFFSPRING = list(bp=bp,cc=cc,mf=mf)

termGene = function( go_id  )
{
        offspringTerm = do.call(c, sapply( GOOFFSPRING , function(x) { as.character(x[,1])[ as.character(x[,2])== go_id ] } ) )
        offspringTerm = c( go_id, offspringTerm )
        term_genes = as.character( annotGO[as.character(annotGO$go_id) %in% offspringTerm, 1])
        unique(term_genes)
}

##############################################################################################
##############################################################################################

setwd('/hpc/users/zhus02/fangg03a/sjzhu/Projects/myISMARA/GSEA/GO')
tag      = unique( paste( forGO[,4] , forGO[,6] , forGO[,8] , sep='\t' ) )
splitTag = strsplit(tag,'\t')
goid     = sapply( splitTag , function(x) x[1]  )
ontology = sapply( splitTag , function(x) x[2]  )
term     = sapply( splitTag , function(x) x[3]  )

geneset = c()
for(  i in 1:length(goid) )
{
  cat(i,'\n')
  geneset[i] = paste( termGene(goid[i]) , collapse=',' )
}

GO.gs = data.frame(goid,ontology,term,geneset)
save(GO.gs,file='GO.gs.rda')

##############################################################################################
##############################################################################################
setwd('/hpc/users/zhus02/fangg03a/sjzhu/Projects/myISMARA/GSEA/GO')
load('GO.gs.rda')
net = geneSet2sparseMatrix( term=GO.gs$goid , geneset=GO.gs$geneset )
annot = GO.gs[,-4]
annot$totalGenes = as.integer(colSums(net))
GO = list(net=net,annot=annot)
save(GO,file='GO.rda')

object.size(GO.gs)/1e6
object.size(GO)/1e6




