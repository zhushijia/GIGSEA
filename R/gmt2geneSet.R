#' gmt2geneSet
#'
#' gmt2geneSet transforms a gmt format file into geneSets.
#'
#' @param gmt a vector of character values. Each item is a list of words comprising a term and its corresponding gene set, which are separated by tab.
#' @param termCol an integer value indicating in each item of gmt, which word is the term , by default, 1.
#' @param nonGeneCol an integer value indicating in each item of gmt, which words are not the gene set, by default, termCol.
#' @param singleValue a numeric value, which assigns the same value to all genes in a given gene set. This is useful when combining together the up-regulated gene sets (regularly, singleValue=1) and the down-regulated gene sets (regularly, singleValue=-1)
#'
#' @return a data frame, comprising three vectors: term (like pathway names), geneset (a gene symbol list separate by comma), and value (either discrete or continuous separated by comma)
#' @export
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @seealso \code{\link{geneSet2Net}}; \code{\link{geneSet2sparseMatrix}};
#'
gmt2geneSet <- function( gmt , termCol=1 , nonGeneCol=termCol , singleValue=NULL )
{
  split_gmt = strsplit( as.character(gmt) , '\t' )
  term = sapply( split_gmt , function(x) paste( x[termCol],collapse=',' ) )
  geneset = sapply( split_gmt , function(x) paste( x[-nonGeneCol],collapse=',' ) )
  if( is.null(singleValue) )
  {
    gff = data.frame(term , geneset , stringsAsFactors = FALSE)
  } else {
    gene_num <- sapply( strsplit(geneset,',') , length )
    value <- sapply( gene_num , function(times) paste( rep(singleValue,times),collapse=',') )
    gff = data.frame(term , geneset , value , stringsAsFactors = FALSE)
  }
  gff
}
