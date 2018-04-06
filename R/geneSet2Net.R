#' geneSet2Net
#'
#' geneSet2Net transforms gene sets to a matrix, which represents the connectivity between terms and genes.
#'
#' @param term a vector of character values incidating the names of gene sets, like pathway names and miRNA names.
#' @param geneset a vector of character values, where each value is a gene list separated by 'sep'.
#' @param value a vector of numeric values indicating the connectivity of between terms and genes. It could take either discrete values (0 and 1) or continuous vlaues.
#' @param sep a character which separates the genes in the geneset.
#'
#' @return
#' a matrix of numeric values where the column corresponds to the term and the row corresponds to the geneset.
#'
#' @export
#'
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#'
#' @seealso \code{\link{geneSet2sparseMatrix}}; \code{\link{gmt2geneSet}};
#'
geneSet2Net <- function( term , geneset , value=NULL , sep=',' )
{

  split_geneset = strsplit( as.character(geneset) , sep )
  names(split_geneset) = term
  all_genes = unique( do.call(c,split_geneset) )

  if( is.null(value) )
  {
    net = sapply( 1:length(term) , function(i) {
      neti = rep(0,length(all_genes))
      neti[ all_genes %in% split_geneset[[i]] ] = 1
      neti
    } )
  } else {
    split_value = lapply( strsplit( as.character(value) , sep ) ,  function(x) as.numeric(as.character(x)) )
    net = sapply( 1:length(term) , function(i) {
      neti = rep(0,length(all_genes))
      index = match( split_geneset[[i]] , all_genes )
      neti[ index  ] = split_value[[i]]
      neti
    } )
  }

  colnames(net) = term
  rownames(net) = all_genes
  net = net[ order(all_genes) , ]
  net

}

