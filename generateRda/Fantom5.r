

setwd('E:/SJZHU/Projects/Heterogeneity/ISMARA_targetScan_refined_pipeline/data/NewData')
motif_count = read.table('hg19_f5_sitecount_matrix',sep='\t',header=T)
promoter_gene = readLines('hg19_f5_prom_annotation')[-1]
pg_promoter = sapply( strsplit(promoter_gene,'\t') , function(x) x[6] )
pg_gene = sapply( strsplit(promoter_gene,'\t') , function(x) {
		strsplit( as.character(x[11]),'[|]')[[1]][2]
		})
promoter_gene = data.frame(promoter=pg_promoter , gene=pg_gene)
promoter_gene = promoter_gene[ match( as.character(motif_count$promoter) , as.character(promoter_gene$promoter) ) ,   ]
all( as.character(motif_count$promoter) == as.character(promoter_gene$promoter)	 )

tf_count = motif_count[,-c(1,502:604)]
group = as.character(promoter_gene$gene)
tf_avgCount = apply( tf_count ,  2 , function(x) { tapply(x,group,mean) } )
# get mean, not all promoter of one gene have the same tf binding, so that tf on some promoters of one gene are 0



# save the data
term = colnames(tf_count)
geneset = apply( tf_avgCount , 2 , function(x) paste(rownames(tf_avgCount)[x>0],collapse=',') )
value = apply( tf_avgCount , 2 , function(x) paste(x[x>0],collapse=',') )

Fantom5.TF.gs = data.frame( term=term , geneset=geneset , value=value )
Fantom5.TF.net = geneSet2sparseMatrix( term=Fantom5.TF.gs$term , geneset=Fantom5.TF.gs$geneset , value=Fantom5.TF.gs$value )

TF = colnames(Fantom5.TF.net)
totalGenes = apply( Fantom5.TF.net , 2 , function(x) sum(x>0,na.rm=T) )
annot = data.frame( TF , totalGenes )

Fantom5.TF = list(net=Fantom5.TF.net,annot=annot)

setwd("C:/Users/zhus02/Dropbox/mypaper/GIGSEA/codes/GIGSEA")
save(Fantom5.TF.gs,file='data/Fantom5.TF.gs.rda')
save(Fantom5.TF,file='data/Fantom5.TF.rda')
