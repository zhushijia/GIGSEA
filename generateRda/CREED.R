gmt <- readLines('http://amp.pharm.mssm.edu/CREEDS/download/single_drug_perturbations-v1.0.gmt')
range_up <- grep('-up',gmt)
range_down <- grep('-dn',gmt)
gff_up <- gmt2gff( gmt[range_up] , termCol=c(1,2) , singleValue = 1 )
gff_down <- gmt2gff( gmt[range_down] , termCol=c(1,2) , singleValue = -1 )

term_up <- sapply( gff_up$term , function(x) gsub('-up','',x) )
term_down <- sapply( gff_down$term , function(x) gsub('-dn','',x) )
all(term_up==term_down)

geneset <- sapply( 1:nrow(gff_up) , function(i) paste(gff_up$geneset[i],gff_down$geneset[i],sep=',') )
value <- sapply( 1:nrow(gff_up) , function(i) paste(gff_up$value[i],gff_down$value[i],sep=',') )
net <- geneSet2sparseMatrix( term=term_up , geneset=geneset , value=value )




