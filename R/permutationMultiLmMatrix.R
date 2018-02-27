#' permutationMultiLmMatrix
#'
#' @param fc 
#' @param net 
#' @param weights 
#' @param num 
#' @param step 
#'
#' @return
#' @export
#'
#' @examples
#' 
permutationMultiLmMatrix = function( fc , net , weights=rep(1,nrow(net)) , num=100 , step=1000 )
{
  fc[is.na(fc)] = 0
  weights[is.na(weights)]=0
  net = net[,colSums(net)>0]
  observedTstats = weightedMultiLm( x=net , y=fc, w=weights )[,1]
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
    shuffledTstats =  weightedMultiLm( x=net , y=shuffledFC, w=weights )
    
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
  empiricalPval = empiricalSum/(num*ncol(net))
  lc = locfdr( -1 * sign(observedTstats) * qnorm(empiricalPval) , plot=0 )
  localFdr = lc$fdr
  p0 = lc$fp0[3,3]
  BayesFactor = (1-localFdr)/localFdr *  p0/(1-p0)
  
  pval = data.frame( term , usedGenes , observedTstats , observedPval , empiricalPval , localFdr , BayesFactor )
  rownames(pval) = NULL
  pval = pval[order(pval$BayesFactor,decreasing=T),]
  pval
  
}

