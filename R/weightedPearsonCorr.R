#' weightedPearsonCorr
#'
#' weightedPearsonCorr caculates the weighted Pearson correlation
#'
#' @param x a matrix of numeric values in the size of genes x featureA
#' @param y a matrix of numeric values in the size of genes x featureB
#' @param w a vector of numeric values indicating the weights of genes
#'
#' @return a matrix of numeric values in the size of featureA*featureB, indicating the weighted Pearson correlation coefficients
#' @export
#'
#' @examples
#'
#'
#' # load data
#' data(heart.metaXcan)
#' gene = heart.metaXcan$gene_name
#'
#' # extract the imputed Z-score of gene differential expression, which follows the normal distribution
#' fc <- heart.metaXcan$zscore
#'
#' # use as weights the prediction R^2 and fraction of imputation-used SNPs 
#' usedFrac <- heart.metaXcan$n_snps_used / heart.metaXcan$n_snps_in_cov
#' r2 <- heart.metaXcan$pred_perf_r2
#' weights <- usedFrac*r2
#'
#' # build a new data frame for the following weighted simple linear regression-based enrichment analysis
#' data <- data.frame(gene,fc,weights)
#' head(data)
#'
#' net <- MSigDB.KEGG.Pathway$net
#'
#' # intersect the imputed genes with the gene sets of interest
#' data2 <- orderedIntersect( x = data , by.x = data$gene , by.y = rownames(net)  )
#' net2 <- orderedIntersect( x = net , by.x = rownames(net) , by.y = data$gene  )
#' all( rownames(net2) == as.character(data2$gene) )
#'
#' # calculate the weighted Pearson correlation
#' observedCorr = weightedPearsonCorr( x=net2 , y=data2$fc, w=data2$weights )
#'
#' # calculate the p values of the weighted Pearson correlation
#' observedPval = matrixPval( observedCorr, df=sum(weights>0,na.rm=TRUE)-2 )
#'
#' res = data.frame( observedCorr , observedPval )
#' head(res)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @seealso \code{\link{orderedIntersect}}; \code{\link{matrixPval}};
#'
weightedPearsonCorr <- function( x, y, w=rep(1,nrow(x))/nrow(x) )
{
  x = as.matrix(x)
  y = as.matrix(y)
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  w[is.na(w)] = 0

  w <- w / sum(w)
  x <- sweep(x, 2, colSums(x * w))
  y <- sweep(y, 2, colSums(y * w))
  t(w*x)%*% y / sqrt( colSums(w * x**2) %*% t(colSums(w * y**2)) )
}


