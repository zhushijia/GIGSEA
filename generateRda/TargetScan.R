library(Matrix)
setwd('E:/SJZHU/Projects/Heterogeneity/ISMARA_targetScan_refined_pipeline/data/HUMAN/miRNA_targetScan')
mirna_count = read.table('Summary_Counts_hsa.txt',sep='\t',header=T)
mirna_count = subset( mirna_count , Aggregate.PCT!="NULL" & Total.num.conserved.sites>0 , select=c(miRNA.family,Gene.Symbol,Transcript.ID,Aggregate.PCT) )
mirna_count$Aggregate.PCT = as.numeric(as.character(mirna_count$Aggregate.PCT))
colnames(mirna_count) = c('family','gene','transcript','pct')
length( unique( subset( mirna_count, family=='CUUUGGU')$gene)) # mir9 targets: both 5 and sum(6:8) have the same mir9 targets (1237)
length(unique(mirna_count$family)) # 87

mirna_family = read.table('miR_Family_Info.txt',sep='\t',header=T)
mirna_family = subset( mirna_family , grepl('hsa',MiRBase.ID) & Family.Conservation.>1 )
mirna_family = tapply( as.character(mirna_family$MiRBase.ID) , as.character(mirna_family$Seed.m8) , function(x) paste(x,collapse=" "))
mirna_family = data.frame( family = names(mirna_family) , mirna=mirna_family )
rownames(mirna_family) = NULL

group = paste( as.character(mirna_count$family) , as.character(mirna_count$gene) , sep = ',' )
avg_pct = tapply( mirna_count$pct , group , mean )
# because the duplicates of genes are for transcripts, so getting mean makes sense
# in addition, almost all transcript of one gene have the same value, so,
mirna_gene = t( sapply( strsplit(names(avg_pct),',') , function(x){x} ) )

mirna_avgCount = data.frame( family = mirna_gene[,1] , gene =mirna_gene[,2] , avg_pct = avg_pct )
mirna_avgCount = merge( mirna_family , mirna_avgCount , by='family' )


# transformation and save
TargetScan.miRNA.gs = dataFrame2geneSet( term=mirna_avgCount$mirna , gene=mirna_avgCount$gene , value=mirna_avgCount$avg_pct)
TargetScan.miRNA.net = geneSet2sparseMatrix( term=TargetScan.miRNA.gs$term , geneset=TargetScan.miRNA.gs$geneset , value = TargetScan.miRNA.gs$value )

annot=mirna_family[ match( colnames(TargetScan.miRNA.net) , as.character(mirna_family[,2]) ) ,  ]
annot = annot[,c(2,1)]
annot$totalGenes = apply( TargetScan.miRNA.net, 2, function(x) sum(x!=0,na.rm=T) )
TargetScan.miRNA = list( net=TargetScan.miRNA.net , annot=annot )

setwd("C:/Users/zhus02/Dropbox/mypaper/GIGSEA/codes/GIGSEA")
save(TargetScan.miRNA.gs,file='data/TargetScan.miRNA.gs.rda')
save(TargetScan.miRNA,file='data/TargetScan.miRNA.rda')



