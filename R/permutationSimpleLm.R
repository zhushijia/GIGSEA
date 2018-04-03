#' permutationSimpleLm
#'
#' permutationSimpleLm is a permutation test to calculate the empirical p values for a weighted simple linear regression.
#'
#' @param fc a vector of numeric values representing the gene expression fold change
#' @param net a matrix of numeric values in the size of gene number x gene set number, representing the connectivity betweeen genes and gene sets
#' @param weights a vector of numeric values representing the weights of permuted genes
#' @param num a vector of integer values representing the number of permutations
#'
#' @return a data frame comprising the following columns:
#' \itemize{
#' \item {term} a vector of character values incidating the name of gene set.
#' \item {usedGenes} a vector of numeric values indicating the number of genes used in the model.
#' \item {Estimate} a vector of numeric values indicating the regression coefficients.
#' \item {Std..Error} a vector of numeric values indicating the standard errors of regression coefficients.
#' \item {t.value} a vector of numeric values indicating the t-statistics of regression coefficients.
#' \item {observedPval} a vector of numeric values [0,1] indicating the p values from weighted simple regression model.
#' \item {empiricalPval} a vector of numeric values [0,1] indicating the empirical p values from the permutation test.
#' }
#'
#'
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
#' # extract the imputed Z-score of differential gene expression, which follows the normal distribution
#' fc <- heart.metaXcan$zscore
#'
#' # use as weights the prediction R^2 and the fraction of imputation-used SNPs 
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
#' # the SGSEA.res1 uses the  weighted simple linear regression model, while SGSEA.res2 used the weighted Pearson correlation. The latter one takes substantially less time.
#' system.time( SGSEA.res1 <- permutationSimpleLm( fc=data2$fc , net=net2 , weights=data2$weights , num=100 ) )
#' system.time( SGSEA.res2 <- permutationSimpleLmMatrix( fc=data2$fc , net=net2 , weights=data2$weights , num=100 ) )
#' head(SGSEA.res2)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{orderedIntersect}}; \code{\link{permutationSimpleLmMatrix}};
#'
permutationSimpleLm = function( fc , net , weights=rep(1,nrow(net)) , num=100 )
{
  shuffledPval <- lapply( 1:num, function(i) {
    reportProgress(i,num,10)
    shuffledFC = sample(fc,length(fc))
    separateLm( shuffledFC , net , weights )
  })

  shuffledPval = do.call( cbind , shuffledPval )
  observedPval = separateLm( fc , net , weights )

  empiricalPval = c()
  for(i in 1:nrow(shuffledPval))
  {
    empiricalPval[i] = mean( shuffledPval[i,]<observedPval[i] )
  }

  term = colnames(net)
  usedGenes = apply(as.matrix(net), 2, function(x) sum(x!=0,na.rm=T) )
  pval = data.frame( term , usedGenes , observedPval , empiricalPval )
  rownames(pval) = NULL
  pval = pval[order(pval$empiricalPval),]
  pval
}


separateLm = function( y , net , weights )
{
  apply( net , 2 , function(x) {
    z = lm( y~x , weights=weights )
    summary(z)$coef[-1,4]
  })
}

