#' empiricalPval
#'
#' empiricalPval estimates the single-sided empirical pvalues of observed 
#' statistics based on the shuffled statistics 
#'
#' @param observedStats a vector of observed statistics
#' @param shuffledStats a vector of shuffled statistics
#'
#' @return a vector of p values taking values in [0,1]
#' @export
#'
#' @examples
#'
#' observedStats = qnorm( seq(0,1,0.1) )
#' shuffledStats = rnorm(100000)
#' # estimate the single-sided empirical pvalue
#' emp = getEmpiricalNum(observedStats,shuffledStats) 
#' cat(round(emp,2),'\n') 
#' 
empiricalPval = function( observedStats , shuffledStats , step=1 )
{
  
  empiricalNum = rep(0,length(observedStats))
  ord = order(observedStats)
  sortedShuffledStats = sort(shuffledStats)
  L = length(sortedShuffledStats)
  
  k = 1
  for( i in seq_len(length(observedStats)) )
  {
    j = ord[i]
    while( sortedShuffledStats[k]<observedStats[j] & k<L ) 
    { 
      k <- k+step 
      if (k>L)
      {
        k <- L
        break
      }
    }
    while( k>0 & sortedShuffledStats[k]>=observedStats[j]  ) 
    { 
      k <- k-1 
      if(k<=0) 
      {
        k <- 1
        break
      }
    }
    if(observedStats[j]>0) 
    {
      empiricalNum[j] <- L - k
    } else {
      empiricalNum[j] <- k
    }
  }
  
  empiricalNum/L
  
}

#system.time(num11<-empiricalPval( observedCorr , shuffledCorr , step=1000 ))

