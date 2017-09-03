buildGeneset <- function( gmt )
{
  split_gmt = strsplit( as.character(gmt) , '\t' )
  term = sapply( split_gmt , function(x) x[1] )
  link = sapply( split_gmt , function(x) x[2] )
  geneset = sapply( split_gmt , function(x) paste( x[-c(1,2)],collapse=',' ) )

  gff = data.frame(term,link,geneset)
  #gff = gff[ order(as.character(gff$term)), ]
  gff
}




