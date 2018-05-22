#' geneSet2Net
#'
#' geneSet2Net transforms gene sets to a matrix, which represents the 
#' connectivity between terms and genes.
#'
#' @param term a vector of character values incidating the names of gene sets, 
#' like pathway names and miRNA names.
#' @param geneset a vector of character values, where each value is a gene 
#' list separated by 'sep'.
#' @param value a vector of numeric values indicating the connectivity of 
#' between terms and genes. It could take either discrete values (0 and 1) 
#' or continuous vlaues.
#' @param sep a character which separates the genes in the geneset.
#'
#' @return
#' a matrix of numeric values where the column corresponds to the term and the 
#' row corresponds to the geneset.
#'
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
#' # combine up and down-regulated gene sets, and use 1 and -1 to indicate 
#' # their direction 
#' # extract the drug names
#' term_up<-vapply( gff_up$term, function(x) gsub('-up','',x), character(1) )
#' term_down<-vapply( gff_down$term, function(x) gsub('-dn','',x), character(1))
#' all(term_up==term_down)
#'
#' # combine the up-regulated and down-regulated gene names for each 
#' # drug perturbation
#' geneset <- vapply(1:nrow(gff_up),function(i) paste(gff_up$geneset[i],
#' gff_down$geneset[i],sep=','), character(1) )
#'
#' # use 1 and -1 to indicate the direction of up and down-regulated genes
#' value <- vapply( 1:nrow(gff_up) , function(i) paste(gff_up$value[i],
#' gff_down$value[i],sep=',') , character(1) )
#'
#' # transform the gene set into matrix, where the row represents the gene, 
#' # the column represents the drug perturbation, and each entry takes values 
#' # of 1 and -1
#' net1 <- geneSet2Net( term=term_up , geneset=geneset , value=value )
#' # transform the gene set into sparse matrix, where the row represents the 
#' # gene, the column represents the drug perturbation, and each entry takes 
#' # values of 1 and -1
#' net2 <- geneSet2sparseMatrix( term=term_up , geneset=geneset , value=value )
#' tail(net1[,1:30])
#' tail(net2[,1:30])
#' # the size of sparse matrix is much smaller than the matrix
#' format( object.size(net1), units = "auto")
#' format( object.size(net2), units = "auto")
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
    net = vapply( seq_along(term) , function(i) {
      neti = rep(0,length(all_genes))
      neti[ all_genes %in% split_geneset[[i]] ] = 1
      neti
    } , numeric(length(all_genes)) )
  } else {
    sv_tmp = strsplit( as.character(value) , sep )
    split_value = lapply( sv_tmp ,  function(x) as.numeric(as.character(x)) )
    net = vapply( seq_along(term) , function(i) {
      neti = rep(0,length(all_genes))
      index = match( split_geneset[[i]] , all_genes )
      neti[ index  ] = split_value[[i]]
      neti
    } , numeric(length(all_genes)) )
  }

  colnames(net) = term
  rownames(net) = all_genes
  net = net[ order(all_genes) , ]
  net

}

