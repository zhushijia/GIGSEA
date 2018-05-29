#' dataframe2geneSet
#'
#' dataframe2geneSet transforms a data frame (1term-1gene) into geneSets 
#' (1term-Ngenes).
#'
#' @param term a character value incidating the name of the column for the 
#' gene set (terms).
#' @param gene a character value incidating the name of the column for the 
#' genes.
#' @param value a vector of numeric values indicating the connectivity of 
#' between terms and genes. It could take either discrete values (0 and 1) 
#' or continuous values.
#'
#' @return a data frame, comprising three vectors: term (like pathway names), 
#' geneset (a gene symbol list separate by comma), and value (either discrete 
#' or continuous separated by comma)
#' 
#' @export
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#'
#' @seealso \code{\link{geneSet2Net}}; \code{\link{geneSet2sparseMatrix}};
#'
dataframe2geneSet <- function( term , gene , value=NULL )
{
    if( is.null(value) )
    {
        geneset <- split( gene , term )
        geneset <- vapply( geneset, function(x) 
            paste(x,collapse=','), character(1) )
        res <- data.frame( term=names(geneset) , geneset=geneset )
    } else {
        gene_value <- paste( gene , value )
        gene_value_set <- split( gene_value , term )
        geneset <- vapply( gene_value_set , function(x) {
            genes <- vapply( strsplit(x,' '), function(gv) gv[1], character(1)) 
            paste(genes,collapse=',')
            } , character(length(gene_value_set)) )
        valueset <- vapply( gene_value_set , function(x) {
            values <- vapply( strsplit(x,' '), function(gv) gv[2], character(1))
            paste( values , collapse=',' )
            } , character(length(gene_value_set)) )
        res <- data.frame( term=names(geneset) , geneset=geneset , 
                        valueset=valueset )
        }
    rownames(res) <- NULL
    res
  
}
