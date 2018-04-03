#' permutationMultipleLm
#'
#' permutationMultipleLm is a permutation test to calculate the empirical p values for a weighted multiple linear regression.
#'
#' @param fc a vector of numeric values representing gene expression fold change
#' @param net a matrix of numeric values in the size of gene number x gene set number, representing the connectivity between genes and gene sets
#' @param weights a vector of numeric values representing the weights of permuated genes
#' @param num an integer value representing the number of permutations
#'
#' @return a data frame comprising the following columns:
#' \itemize{
#' \item {term} a vector of character incidating the names of gene sets.
#' \item {usedGenes} a vector of numeric values indicating the number of genes used in the model.
#' \item {Estimate} a vector of numeric values indicating the regression coefficients.
#' \item {Std..Error} a vector of numeric values indicating the standard errors of regression coefficients.
#' \item {t.value} a vector of numeric values indicating the t-statistics of regression coefficients.
#' \item {observedPval} a vector of numeric values [0,1] indicating the p values from the multiple weighted regression model.
#' \item {empiricalPval} a vector of numeric values [0,1] indicating the empirical p values from the permutation test.
#' }
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
#' # extract the imputed Z-score of the differential gene expression, which follows the normal distribution
#' fc <- heart.metaXcan$zscore
#'
#' # use as weights the prediction R^2 and the fraction of imputation-used SNPs 
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
#' # intersect the imputed genes with the gene sets of interest
#' data2 <- orderedIntersect( x = data , by.x = data$gene , by.y = rownames(net)  )
#' net2 <- orderedIntersect( x = net , by.x = rownames(net) , by.y = data$gene  )
#' all( rownames(net2) == as.character(data2$gene) )
#'
#' # run the MGSEA
#' MGSEA.res <- permutationMultipleLm( fc=data2$fc , net=net2 , weights=data2$weights , num=1000 )
#' head(MGSEA.res)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{orderedIntersect}};
#'
permutationMultipleLm = function( fc , net , weights=rep(1,nrow(net)) , num=100 )
{
  net = as.matrix(net) # to transform the sparseMatrix to matrix
  shuffle <- lapply( 1:num, function(i) {
    reportProgress(i,num,10)
    #randomfc = rnorm(nrow(net))
    shuffledFC = sample(fc,length(fc))
    shuffledLM = lm( shuffledFC~net , weights=weights )
    p = summary(shuffledLM)$coef[-1,4]
    f = summary(shuffledLM)$fstatistic[1]
    list( p=p , f=f )
  } )

  shuffledPval = do.call( cbind , lapply(shuffle,function(x) x$p ) )
  shuffledF    = do.call(   c   , lapply(shuffle,function(x) x$f ) )

  trueLM = lm( fc~net , weights=weights )
  observedCoef = summary(trueLM)$coef[-1,-4]
  observedPval = summary(trueLM)$coef[-1,4]

  observedF = summary(trueLM)$fstatistic[1]

  empiricalPval = c()
  for(i in 1:nrow(shuffledPval))
  {
    empiricalPval[i] = mean( shuffledPval[i,]<observedPval[i] )
  }

  nonNaRange = which( !is.na(trueLM$coefficients[-1]) )
  term = colnames(net)[nonNaRange]
  usedGenes = apply(as.matrix(net), 2, function(x) sum(x!=0,na.rm=T) )[nonNaRange]
  # usedWeightedGenes = colSums( net*weights )
  pval = data.frame( term , usedGenes , observedCoef , observedPval , empiricalPval )
  rownames(pval) = NULL
  pval = pval[order(pval$empiricalPval),]

  fval = c( f=observedF , f.pval=mean( shuffledF<observedF ) )

  #list( fval=fval , pval=pval )

  pval

}
