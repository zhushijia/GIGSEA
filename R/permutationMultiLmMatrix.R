#' permutationMultiLmMatrix
#'
#' permutationMultiLmMatrix is a permutation test to calculate the empirical p values for weighted multiple linear regression
#'
#' @param fc a vector of numeric values representing the gene expression fold change
#' @param net a matrix of numeric values in the size of sample*feature representing the gene sets
#' @param weights a vector of numeric values representing the weights of samples
#' @param num an integer value representing the number of permutations
#' @param step an integer value representing the number of permutations in each step 
#'
#' @return a data frame comprising following columns:
#' \itemize{
#' \item {term} a vector of character incidating the name of gene set.
#' \item {usedGenes} a vector of numeric values indicating the number of gene used in the model.
#' \item {observedTstats} a vector of numeric values indicating the observed t-statistics for the weighted regression coefficients.
#' \item {observedPval} a vector of numeric values [0,1] indicating the observed p values of the weighted regression coefficients.
#' \item {empiricalPval} a vector of numeric values [0,1] indicating the permutation test-resulting p values of the weighted regression coefficients.
#' \item {empiricalPval} a vector of numeric values indicating the Bayes Factor for multiple test correction.
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
#' # the MGSEA.res1 uses the single weighted linear regression model, while MGSEA.res2 used the weighted Pearson correlation. The latter one takes substantially less time.
#' system.time( MGSEA.res1 <- permutationMultiLm( fc=data2$fc , net=net2 , weights=data2$weights , num=100 ) )
#' system.time( MGSEA.res2 <- permutationMultiLmMatrix( fc=data2$fc , net=net2 , weights=data2$weights , num=100 ) )
#' head(MGSEA.res2)
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{orderedIntersect}}; \code{\link{permutationMultiLm}};
#'
permutationMultiLmMatrix = function( fc , net , weights=rep(1,nrow(net)) , num=100 , step=1000 )
{
  fc[is.na(fc)] = 0
  weights[is.na(weights)]=0
  net = as.matrix(net)
  net[is.na(net)] = 0
  net = net[,colSums(net)>0]
  
  x = net
  w = weights
  y = fc
  
  X = cbind(1,x)
  W = diag(w)
  A0 = t(X) %*% W
  A = ginv(A0 %*% X)
  B = A0 %*% y
  coefs = A %*% B
  predicts = X %*% coefs
  residuals = y - predicts
  delta = colSums( w * residuals^2 )/( sum(w>0,na.rm=T)-qr(X)$rank )
  se = sapply( delta , function(d) sqrt( d * diag(A) ) )
  t = coefs/se
  observedTstats = as.matrix(t[-1,1])
  observedPval = 2 * pt(abs(observedTstats),df=sum(weights>0,na.rm=T)-2,lower.tail=FALSE)
  
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
    
    B = A0 %*% shuffledFC
    coefs = A %*% B
    predicts = X %*% coefs
    residuals = shuffledFC - predicts
    delta = colSums( w * residuals^2 )/( sum(w>0,na.rm=T)-qr(X)$rank )
    se = sapply( delta , function(d) sqrt( d * diag(A) ) )
    t = coefs/se
    shuffledTstats = as.matrix(t[-1,])
    
    for(j in 1:nrow(shuffledTstats))
    {
      if(observedTstats[j]>0) 
      {
        empiricalSum[j] = empiricalSum[j] + sum( shuffledTstats > observedTstats[j] )
      } else {
        empiricalSum[j] = empiricalSum[j] + sum( shuffledTstats < observedTstats[j] )
      }
    }
    
    for(j in cs_steps[i]:cs_steps[i+1]  ) reportProgress(j,num,10)
  }
  
  term = colnames(net)
  usedGenes = apply(as.matrix(net), 2, function(x) sum(x!=0,na.rm=T) )
  empiricalPval = (empiricalSum+0.1)/(num*ncol(net))
  lc = locfdr( sign(observedTstats) * qnorm(empiricalPval) , plot=0 )
  #lc = locfdr( observedTstats , plot=0 )
  localFdr = lc$fdr
  p0 = lc$fp0[3,3] - 1.96*lc$fp0[4,3]
  BayesFactor = (1-localFdr)/localFdr *  p0/(1-p0)
  
  #pval = data.frame( term , usedGenes , observedTstats , observedPval , empiricalPval , localFdr , BayesFactor )
  pval = data.frame( term , usedGenes , observedTstats , observedPval , empiricalPval , BayesFactor )
  rownames(pval) = NULL
  pval = pval[order(pval$BayesFactor,decreasing=T),]
  #pval = pval[order(pval$empiricalPval),]
  pval
  
}

