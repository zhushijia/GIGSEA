#' permutationSingleLmMatrix
#'
#' permutationSingleLmMatrix is a permutation test to calculate the empirical p values for weighted Pearson correlation.
#'
#' @param fc a vector of numeric values representing the gene expression fold change
#' @param net a matrix of numeric values in the size of sample*feature representing the gene sets
#' @param weights a vector of numeric values representing the weights of samples
#' @param num a vector of integer values representing the number of permutations
#'
#' @return a data frame comprising following columns:
#' \itemize{
#' \item {term} a vector of character incidating the name of gene set.
#' \item {usedGenes} a vector of numeric values indicating the number of gene used in the model.
#' \item {observedCorr} a vector of numeric values indicating the observed weighted Pearson correlation coefficients.
#' \item {observedPval} a vector of numeric values [0,1] indicating the observed p values of the weighted Pearson correlation coefficients.
#' \item {empiricalPval} a vector of numeric values [0,1] indicating the permutation test-resulting p values of the weighted Pearson correlation coefficients.
#' }
#'
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
#' # the SGSEA.res1 uses the single weighted linear regression model, while SGSEA.res2 used the weighted Pearson correlation. The latter one takes substantially less time.
#' system.time( SGSEA.res1 <- permutationSingleLm( fc=data2$fc , net=net2 , weights=data2$weights , num=100 ) )
#' system.time( SGSEA.res2 <- permutationSingleLmMatrix( fc=data2$fc , net=net2 , weights=data2$weights , num=100 ) )
#' head(SGSEA.res2)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{orderedIntersect}}; \code{\link{permutationSingleLm}};
#'
permutationSingleLmMatrix = function( fc , net , weights=rep(1,nrow(net)) , num=100 , step=1000 )
{

  observedCorr = weightedPearsonCorr( a=fc, b=net , w=weights )[1,]
  observedPval = matrixPval( observedCorr, df=sum(weights>0,na.rm=T)-2 )
  #observedPval2 = separateLm( fc , net , weights )
  #max(abs(observedPval-observedPval2))

  empiricalSum = rep(0,ncol(net))
  if( num%%step==0 )
  {
    steps = rep(step,num%/%step)
  } else {
    steps = c( rep(step,num%/%step) , num%%step )
  }
  cs_steps = c( 1, cumsum(steps) )

  for( i in 1:length(steps) )
  {
    stepi = steps[i]
    shuffledFC = sapply(1:stepi,function(s) sample(fc) )
    shuffledCorr =  weightedPearsonCorr( a=shuffledFC, b=net , w=weights )

    for(j in 1:ncol(shuffledCorr))
    {
      empiricalSum[j] = empiricalSum[j] + sum( abs(shuffledCorr[,j]) > abs(observedCorr[j]) )
    }

    for(j in cs_steps[i]:cs_steps[i+1]  ) reportProgress(j,num,10)
  }

  term = colnames(net)
  usedGenes = apply(as.matrix(net), 2, function(x) sum(x!=0,na.rm=T) )
  empiricalPval = empiricalSum/num
  pval = data.frame( term , usedGenes , observedCorr , observedPval , empiricalPval )
  rownames(pval) = NULL
  pval = pval[order(pval$empiricalPval),]
  pval

}


