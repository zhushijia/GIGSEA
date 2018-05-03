#' permutationSimpleLmMatrix
#'
#' permutationSimpleLmMatrix is a permutation test to calculate the empirical p values for the weighted simple linear regression model based on the weighted Pearson correlation.
#'
#' @param fc a vector of numeric values representing the gene expression fold change
#' @param net a matrix of numeric values in the size of gene number x gene set number, representing the connectivity betwen genes and gene sets
#' @param weights a vector of numeric values representing the weights of permuted genes
#' @param num an integer value representing the number of permutations
#' @param step an integer value representing the number of permutations in each step 
#'
#' @return a data frame comprising following columns:
#' \itemize{
#' \item {term} a vector of character values incidating the name of gene set.
#' \item {usedGenes} a vector of numeric values indicating the number of gene used in the model.
#' \item {observedCorr} a vector of numeric values indicating the observed weighted Pearson correlation coefficients.
#' \item {empiricalPval} a vector of numeric values [0,1] indicating the permutation-based empirical p values. 
#' \item {BayesFactor} a vector of numeric values indicating the Bayes Factor for the multiple test correction.
#' }
#'
#' @export
#'
#' @importFrom stats lm pt qnorm
#' @importFrom utils read.table write.table
#' @importFrom Matrix sparseMatrix
#' @importFrom locfdr locfdr
#' @import MASS
#' 
#' @examples
#'
#' # load data
#' data(heart.metaXcan)
#' gene = heart.metaXcan$gene_name
#'
#' # extract the imputed Z-score of gene differential expression, which follows the normal distribution
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
#' # intersect the permuted genes with the gene sets of interest
#' data2 <- orderedIntersect( x = data , by.x = data$gene , by.y = rownames(net)  )
#' net2 <- orderedIntersect( x = net , by.x = rownames(net) , by.y = data$gene  )
#' all( rownames(net2) == as.character(data2$gene) )
#'
#' # the SGSEA.res1 uses the weighted simple linear regression model, 
#' # while SGSEA.res2 used the weighted Pearson correlation. The latter one takes substantially less time.
#' #system.time(SGSEA.res1<-permutationSimpleLm(fc=data2$fc, net=net2, weights=data2$weights, num=1000))
#' system.time(SGSEA.res2<-permutationSimpleLmMatrix(fc=data2$fc, net=net2, weights=data2$weights, num=1000))
#' head(SGSEA.res2)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @seealso \code{\link{orderedIntersect}}; \code{\link{permutationSimpleLm}};
#'
permutationSimpleLmMatrix = function( fc , net , weights=rep(1,nrow(net)) , num=100 , step=1000 )
{
  fc[is.na(fc)] = 0
  weights[is.na(weights)]=0
  net = as.matrix(net)
  net = net[,colSums(net)>0]
  observedCorr = weightedPearsonCorr( x=net, y=fc, w=weights )[,1]
  observedPval = matrixPval( observedCorr, df=sum(weights>0, na.rm=TRUE)-2 )
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
    shuffledCorr =  weightedPearsonCorr( x=net, y=shuffledFC, w=weights )
    
    for(j in 1:nrow(shuffledCorr))
    {
      if(observedCorr[j]>0) 
      {
        empiricalSum[j] = empiricalSum[j] + sum( shuffledCorr > observedCorr[j] )
      } else {
        empiricalSum[j] = empiricalSum[j] + sum( shuffledCorr < observedCorr[j] )
      }
    }
    
    for(j in cs_steps[i]:cs_steps[i+1]  ) reportProgress(j,num,10)
  }
  
  term = colnames(net)
  usedGenes = apply(as.matrix(net), 2, function(x) sum(x!=0,na.rm=TRUE) )
  empiricalPval = (empiricalSum+0.1)/(num*ncol(net))
  #lc = locfdr( -1 * sign(observedCorr)*qnorm(empiricalPval) , plot=0 )
  #lc = locfdr( observedCorr , plot=0 )
  zscore = qnorm(empiricalPval)
  lc = locfdr( c(zscore,-zscore) , plot=0 )
  localFdr = lc$fdr[1:length(zscore)]
  
  p0 = lc$fp0[3,3]
  # p0 = lc$fp0[3,3] - 1.96*lc$fp0[4,3]
  
  BayesFactor = (1-localFdr)/localFdr *  p0/abs(1-p0)
  
  #pval = data.frame( term , usedGenes , observedCorr , observedPval , empiricalPval , localFdr , BayesFactor )
  pval = data.frame( term , usedGenes , observedCorr , empiricalPval , BayesFactor )
  rownames(pval) = NULL
  #pval = pval[order(pval$BayesFactor,decreasing=TRUE),]
  pval = pval[order(pval$empiricalPval),]
  pval
  
}

