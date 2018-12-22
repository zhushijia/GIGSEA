#' orderedIntersect
#'
#' orderedIntersect sorts a data frame based on a given collumn and intersects 
#' with another vector.
#'
#' @param x a data frame comprising of numeric values
#' @param by.x a vector of character values. The data frame is sorted based on 
#' by.x
#' @param by.y a vector of character values. After being sorted, the rows of x 
#' are further filtered by intersecting by.x with by.y
#' @param duplicate a character value taking values from "mean", "median", 
#' "max", "max", "min", and "remove". It represents how to deal with the 
#' duplicates in x in terms of by.x. 
#'
#' @return a data frame sorted by "by.x" and intersected with "by.y"
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
#' # normal distribution
#' fc <- heart.metaXcan$zscore
#'
#' # use the prediction R^2 and fraction of imputation-used SNPs as weights
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
#' # intersect the user-provided imputed genes with the gene set of interest
#' data2 <- orderedIntersect( x=data[,c("fc","weights")], 
#' by.x=data$gene, by.y=rownames(net) )
#' net2 <- orderedIntersect( x=net, by.x=rownames(net), by.y=data$gene )
#' all( rownames(net2) == rownames(data2) )
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#'
orderedIntersect <- function(x, by.x, by.y, duplicate='mean')
{
    by.x <- as.character(by.x)
    by.y <- as.character(by.y)
    x <- x[ by.x %in% by.y , ]
    by.x <- by.x[ by.x %in% by.y ]
    
    if( any(table(by.x)>1) )
    {
        x.class = class(x)
      
        if(duplicate=='mean') 
            x <- do.call(cbind, lapply(x, function(xi) tapply(xi,by.x,mean) ) )
        if(duplicate=='median') 
            x <- do.call(cbind, lapply(x, function(xi) tapply(xi,by.x,median)))
        if(duplicate=='max') 
            x <- do.call(cbind, lapply(x, function(xi) tapply(xi,by.x,max) ) )
        if(duplicate=='min') 
            x <- do.call(cbind, lapply(x, function(xi) tapply(xi,by.x,min) ) )
        if(duplicate=='remove') 
            x <- x[ ! by.x %in% names(which(table(by.x)>1)) ,  ]
      
        by.x <- rownames(x)
        if( x.class=='data.frame' ) x <- data.frame(x)
    }
    
    rownames(x) <- by.x
    x[order(rownames(x)), ]
}
