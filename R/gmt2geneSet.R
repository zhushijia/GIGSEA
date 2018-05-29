#' gmt2geneSet
#'
#' gmt2geneSet transforms a gmt format file into geneSets.
#'
#' @param gmt a vector of character values. Each item is a list of words 
#' comprising a term and its corresponding gene set, which are separated by tab.
#' @param termCol an integer value indicating in each item of gmt, which word 
#' is the term , by default, 1.
#' @param nonGeneCol an integer value indicating in each item of gmt, which 
#' words are not the gene set, by default, termCol.
#' @param singleValue a numeric value, which assigns the same value to all 
#' genes in a given gene set. This is useful when combining together the 
#' up-regulated gene sets (regularly, singleValue=1) and the down-regulated 
#' gene sets (regularly, singleValue=-1)
#'
#' @return a data frame, comprising three vectors: term (like pathway names), 
#' geneset (a gene symbol list separate by comma), and value (either discrete 
#' or continuous separated by comma)
#' @export
#'
#' @examples
#'
#' # download the gmt file
#' gmt <- readLines( paste0('http://amp.pharm.mssm.edu/CREEDS/download/',
#' 'single_drug_perturbations-v1.0.gmt') ) 
#'
#' # obtain the index of up-regulated and down-regulated gene sets
#' index_up <- grep('-up',gmt)
#' index_down <- grep('-dn',gmt)
#'
#' # transform the gmt file into gene sets. The gene set is a data frame, 
#' # comprising three vectors: 
#' # term (here is drug), geneset (a gene symbol list separate by comma), 
#' # and value (1 and -1 separate by comma)
#' gff_up <- gmt2geneSet( gmt[index_up], termCol=c(1,2), singleValue = 1 )
#' gff_down <- gmt2geneSet( gmt[index_down], termCol=c(1,2), singleValue = -1 )
#'
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @seealso \code{\link{geneSet2Net}}; \code{\link{geneSet2sparseMatrix}};
#'
gmt2geneSet <- function( gmt, termCol=1, nonGeneCol=termCol, singleValue=NULL )
{
    split_gmt <- strsplit( as.character(gmt) , '\t' )
    term <- vapply( split_gmt , function(x) 
        paste( x[termCol],collapse=',' ) , character(1) )
    geneset <- vapply(split_gmt, function(x) 
        paste(x[-nonGeneCol],collapse=',') , character(1) )
    if( is.null(singleValue) )
    {
        gff <- data.frame(term , geneset , stringsAsFactors = FALSE)
    } else {
        gene_num <- vapply( strsplit(geneset,',') , length , integer(1) )
        value <- vapply( gene_num , function(times) 
            paste( rep(singleValue,times),collapse=',') , character(1) )
        gff <- data.frame(term , geneset , value , stringsAsFactors = FALSE)
    }
    gff
}
