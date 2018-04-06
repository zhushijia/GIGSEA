#' geneSet2sparseMatrix
#'
#' geneSet2sparseMatrix transforms gene sets to a sparse matrix, which represents the connectivity between terms and genes.
#'
#'
#' @param term a vector of character values incidating the names of gene sets, e.g., pathway names and miRNA names.
#' @param geneset a vector of character values, where each value is a gene list separated by 'sep'.
#' @param value a vector of numeric values indicating the connectivity of between terms and genes. It could take either discrete values (0 and 1) or continuous values.
#' @param sep a character which separates the genes in the geneset.
#'
#' @return
#' a sparse matrix where the column corresponds to the term and the row corresponds to the geneset.
#'
#' @export
#'
#' @examples
#'
#' library(GIGSEA)
#'
#' # download the gmt file
#' gmt <- readLines('http://amp.pharm.mssm.edu/CREEDS/download/single_drug_perturbations-v1.0.gmt')
#'
#' # obtain the index of up-regulated and down-regulated gene sets
#' index_up <- grep('-up',gmt)
#' index_down <- grep('-dn',gmt)
#'
#' # transform the gmt file into gene sets. The gene set is a data frame, comprising three vectors: 
#' # term (here is drug), geneset (a gene symbol list separate by comma), 
#' # and value (1 and -1 separate by comma)
#' gff_up <- gmt2geneSet( gmt[index_up] , termCol=c(1,2) , singleValue = 1 )
#' gff_down <- gmt2geneSet( gmt[index_down] , termCol=c(1,2) , singleValue = -1 )
#'
#' # combine up and down-regulated gene sets, and use 1 and -1 to indicate their direction 
#' # extract the drug names
#' term_up <- sapply( gff_up$term , function(x) gsub('-up','',x) )
#' term_down <- sapply( gff_down$term , function(x) gsub('-dn','',x) )
#' all(term_up==term_down)
#'
#' # combine the up-regulated and down-regulated gene names for each drug perturbation
#' geneset<-sapply(1:nrow(gff_up),function(i) paste(gff_up$geneset[i],gff_down$geneset[i],sep=','))
#'
#' # use 1 and -1 to indicate the direction of the up and down-regulated genes, respectively
#' value <- sapply( 1:nrow(gff_up) , function(i) paste(gff_up$value[i],gff_down$value[i],sep=','))
#'
#' # transform the gene set into matrix, where the row represents the gene, 
#' # the column represents the drug perturbation, and each entry takes values of 1 and -1
#' net1 <- geneSet2Net( term=term_up , geneset=geneset , value=value )
#' # transform the gene set into sparse matrix, where the row represents the gene, 
#' # the column represents the drug perturbation, and each entry takes values of 1 and -1
#' net2 <- geneSet2sparseMatrix( term=term_up , geneset=geneset , value=value )
#' tail(net1[,1:30])
#' tail(net2[,1:30])
#' # the size of sparse matrix is much smaller than the matrix
#' format( object.size(net1), units = "auto")
#' format( object.size(net2), units = "auto")
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#'
#' @seealso \code{\link{gmt2geneSet}}; \code{\link{geneSet2Net}};
#'
geneSet2sparseMatrix <- function( term , geneset , value=NULL , sep=',' )
{
  library(Matrix)
  split_geneset = strsplit( as.character(geneset) , sep )
  names(split_geneset) = as.character(term)
  genes = unique( do.call(c,split_geneset) )

  index = lapply( split_geneset , function(x) match( x , genes )  )
  num = sapply(index,length)
  i = do.call(c,index)
  j = rep(1:length(term),num)

  if( is.null(value) )
  {
    net = sparseMatrix( i , j , x=1 , dims=c(length(genes),length(term)) )
  } else {
    split_value = lapply( strsplit( as.character(value) , sep ) ,  function(x) as.numeric(as.character(x)) )
    x = do.call(c,split_value)
    net = sparseMatrix( i , j , x=x , dims=c(length(genes),length(term)) )
  }

  colnames(net) = term
  rownames(net) = genes
  net = net[ order(genes) , ]
  net

}
