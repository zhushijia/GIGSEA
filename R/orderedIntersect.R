#' orderedIntersect
#'
#' orderedIntersect sorts a data frame based on a given collumn and intersects with another vector.
#'
#' @param x a data frame
#' @param by.x a vector of character values. The data frame is sorted based on by.x
#' @param by.y a vector of character values. After being sorted, the rows of x are further filtered by intersecting by.x with by.y
#'
#' @return a data frame sorted by "by.x" and intersected with "by.y"
#' @export
#'
#' @examples
#'
#' library(GIGSEA)
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
#' # intersect the user-provided imputed genes with the gene set of interest
#' data2 <- orderedIntersect( x = data , by.x = data$gene , by.y = rownames(net)  )
#' net2 <- orderedIntersect( x = net , by.x = rownames(net) , by.y = data$gene  )
#' all( rownames(net2) == as.character(data2$gene) )
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
orderedIntersect <- function(x, by.x, by.y)
{
  by.x <- as.character(by.x)
  by.y <- as.character(by.y)
  x <- x[ by.x %in% by.y , ]
  by.x <- by.x[ by.x %in% by.y ]
  x[order(by.x), ]
}
