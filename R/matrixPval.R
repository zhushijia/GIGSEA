#' matrixPval
#'
#' matrixPval calculates the p values for the correlation coefficients based on t-statistics
#'
#' @param r a vector or matrix of correlation coefficients in [-1,+1]
#' @param df the degree of freedom
#'
#' @return a vector or matrix of p values in [0,1]
#' @export
#'
#' @examples
#'
#' library(GIGSEA)
#' library(Matrix)
#' r <- cor(USArrests)
#' df <- nrow(USArrests) - 2
#' pval1 <- matrixPval(r,df)
#'
#' pval2 <- matrix(ncol=ncol(USArrests),nrow=ncol(USArrests),data=0)
#' for(i in 1:ncol(USArrests))
#' {
#'    for(j in 1:ncol(USArrests))
#'    {
#'      pval2[i,j] <- cor.test(USArrests[,i],USArrests[,j])$p.val
#'    }
#' }
#'
#' head(pval1)
#' head(pval2)
#'
#'
matrixPval = function(r,df)
{
  t = sqrt(df)*abs(r) / sqrt(1-r^2)
  2 * pt(t,df,lower.tail=FALSE)
}
