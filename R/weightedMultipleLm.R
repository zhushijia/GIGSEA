#' weightedMultipleLm
#'
#' weightedMultipleLm solves the weighted multiple linear regression model via 
#' matrix operation
#'
#' @param x a matrix of numeric values in the size of genes x featureA
#' @param y a matrix of numeric values in the size of genes x featureB
#' @param w a vector of numeric values indicating the weights of genes
#'
#' @return a matrix of numeric values in the size of featureA*featureB, 
#' indicating the weighted multiple regression coefficients
#' @export
#'
#' @examples
#'
#'
#' # load data
#' data(heart.metaXcan)
#' gene <- heart.metaXcan$gene_name
#'
#' # extract the imputed Z-score of gene differential expression, which follows 
#' # the normal distribution
#' fc <- heart.metaXcan$zscore
#'
#' # use as weights the prediction R^2 and the fraction of imputation-used SNPs 
#' usedFrac <- heart.metaXcan$n_snps_used / heart.metaXcan$n_snps_in_cov
#' r2 <- heart.metaXcan$pred_perf_r2
#' weights <- usedFrac*r2
#'
#' # build a new data frame for the following weighted linear regression-based 
#' # enrichment analysis
#' data <- data.frame(gene,fc,weights)
#' head(data)
#'
#' net <- MSigDB.KEGG.Pathway$net
#'
#' # intersect the permuated genes with the gene sets of interest
#' data2 <- orderedIntersect( x = data , by.x = data$gene , 
#' by.y = rownames(net)  )
#' net2 <- orderedIntersect( x = net , by.x = rownames(net) , 
#' by.y = data$gene  )
#' all( rownames(net2) == as.character(data2$gene) )
#'
#' # perform the weighted multiple linear regression 
#' observedTstats = weightedMultipleLm( x=net2 , y=data2$fc, w=data2$weights )
#'
#' # calculate the p values of the weighted multiple regression coefficients
#' observedPval = 2 * pt(abs(observedTstats), df=sum(weights>0,na.rm=TRUE)-2, 
#' lower.tail=FALSE)
#'
#' res = data.frame( observedTstats , observedPval )
#' head(res)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @seealso \code{\link{orderedIntersect}}; \code{\link{matrixPval}};
#'
weightedMultipleLm <- function( x, y, w=rep(1,nrow(x))/nrow(x) )
{
  #library(MASS)
  
    x <- as.matrix(x)
    y <- as.matrix(y)
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
    w[is.na(w)] <- 0

    X <- cbind(1,x)
    W <- diag(w)
    #A <- solve(t(X) %*% W %*% X)
    A <- ginv(t(X) %*% W %*% X)
    B <- t(X) %*% W %*% y
    coefs <- A %*% B
    predicts <- X %*% coefs
    residuals <- y - predicts
    #delta <- colSums( w * residuals^2 )/( sum(w>0,na.rm=T)-ncol(X) )
    delta <- colSums( w * residuals^2 )/( sum(w>0,na.rm=TRUE)-qr(X)$rank )
  
    se <- vapply( delta , function(d) sqrt( d * diag(A) ) , numeric(nrow(A)) )

    t <- coefs/se
    as.matrix(t[-1,])
  
}


