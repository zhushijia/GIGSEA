dataFrame2geneSet = function( term , gene , value=NULL , sep=',' )
{
	if( is.null(value) )
	{
		geneset = split( gene , term )
		geneset = sapply( geneset , function(x) paste(x,collapse=',') )
		res = data.frame( term=names(geneset) , geneset=geneset )
	} else {
	
		gene_value = paste( gene , value )
		gene_value_set = split( gene_value , term )
		geneset = sapply( gene_value_set , function(x) {
			genes = sapply( strsplit( x , ' ' ) , function(gv) gv[1] )
			paste(genes,collapse=',')
		} )
		valueset = sapply( gene_value_set , function(x) {
			values = sapply( strsplit( x , ' ' ) , function(gv) gv[2] )
			paste(values,collapse=sep)
		} )
		res = data.frame( term=names(geneset) , geneset=geneset , valueset=valueset )
	}

	rownames(res) = NULL
	res
	
}
