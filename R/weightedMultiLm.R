#' weightedMultiLm
#'
#' weightedMultiLm caculates the weighted Pearson correlation
#'
#' @param a a matrix of numeric values in the size of sample*featureA
#' @param b a matrix of numeric values in the size of sample*featureB
#' @param w a vector of numeric values indicating the weights of samples
#'
#' @return a matrix of numeric values in the size of featureA*featureB, indicating the weighted Pearson correlation coefficients
#' @export
#'
#' @examples
#'
#' library(GIGSEA)
#' library(Matrix)
#'
#' # load data
#' data(heart.metaXcan)
#' gene = heart.metaXcan$gene_name
#'
#' # extract the imputed Z-score of gene differential expression, which follows normal distribution
#' fc <- heart.metaXcan$zscore
#'
#' # use the prediction R^2 and fraction of imputation-used SNPs as weights
#' usedFrac <- heart.metaXcan$n_snps_used / heart.metaXcan$n_snps_in_cov
#' r2 <- heart.metaXcan$pred_perf_r2
#' weights <- usedFrac*r2
#'
#' # build a new data frame for the following weighted linear regression-based enrichment analysis
#' data <- data.frame(gene,fc,weights)
#' head(data)
#'
#' net <- MSigDB.KEGG.Pathway$net
#'
#' # do intersection of genes between the user-provided imputed gene expression dataset and the gene sets of interest
#' data2 <- orderedIntersect( x = data , by.x = data$gene , by.y = rownames(net)  )
#' net2 <- orderedIntersect( x = net , by.x = rownames(net) , by.y = data$gene  )
#' all( rownames(net2) == as.character(data2$gene) )
#'
#' # calculate the weighted Pearson correlation
#' wCorr = weightedPearsonCorr( a=data2$fc, b=net2 , w=data2$weights )[1,]
#'
#' # calculate the p values of the weighted Pearson correlation
#' Pval = matrixPval( observedCorr, df=sum(weights>0,na.rm=T)-2 )
#'
#' res = data.frame( wCorr , Pval )
#' head(res)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{orderedIntersect}}; \code{\link{matrixPval}};
#'
weightedMultiLm <- function( x, y, w=rep(1,nrow(x))/nrow(x) )
{
  x = as.matrix(x)
  y = as.matrix(y)
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  w[is.na(w)] = 0

  X = cbind(1,x)
  W = diag(w)
  A = solve(t(X) %*% W %*% X)
  B = t(X) %*% W %*% y
  coefs = A %*% B
  predicts = X %*% coefs
  residuals = y - predicts
  delta = colSums( w * residuals^2 )/( sum(weights>0,na.rm=T)-ncol(X) )
  se = sapply( delta , function(d) sqrt( d * diag(A) ) )

  t = coefs/se
  as.matrix(t[-1,])
  
}


