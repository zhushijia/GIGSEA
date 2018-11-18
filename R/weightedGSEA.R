#' weightedGSEA
#'
#' weightedGSEA performs both SGSEA and MGSEA for a given list of gene sets, 
#' and writes out the results.
#'
#' @param data a data frame comprising comlumns: gene names (characer), 
#' differential gene expression (numeric) and permuated gene weights (numeric 
#' and optional)
#' @param geneCol an integer or a character value indicating the column of gene 
#' name
#' @param fcCol an integer or a character value indicating the column of 
#' differential gene expression
#' @param weightCol an integer or a character value indicating the column of 
#' gene weights
#' @param geneSet a vector of character values indicating the gene sets of 
#' interest. 
#' @param permutationNum an integer value indicating the number of permutation
#' @param outputDir a character value indicating the directory for saving the 
#' results
#' @param MGSEAthres an integer value indicating the thresfold for MGSEA. MGSEA 
#' is performed with no more than "MGSEAthres" gene sets 
#' @param verbose an boolean value indicating whether or not to print output to 
#' the screen 
#' 
#' 
#' @return TRUE
#' @export
#'
#' @examples
#'
#' data(heart.metaXcan)
#' gene <- heart.metaXcan$gene_name
#' fc <- heart.metaXcan$zscore
#' usedFrac <- heart.metaXcan$n_snps_used / heart.metaXcan$n_snps_in_cov
#' r2 <- heart.metaXcan$pred_perf_r2
#' weights <- usedFrac*r2
#' data <- data.frame(gene,fc,weights)
#' # run one-step GIGSEA 
#' # weightedGSEA(data, geneCol='gene', fcCol='fc', weightCol= 'weights', 
#' #  geneSet=c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA") ,
#' #  permutationNum=10000, outputDir="./GIGSEA" )
#' # dir("./GIGSEA")
#' 
weightedGSEA <- function( data , geneCol , fcCol , weightCol=NULL ,
                        geneSet=c("MSigDB.KEGG.Pathway","MSigDB.TF",
                                  "MSigDB.miRNA") ,
                        permutationNum=100 , outputDir=getwd() , 
                        MGSEAthres = NULL , verbose = TRUE )
{
  
    if( !file.exists(outputDir) )
    { 
        if( verbose ) message('creating ' , outputDir )
        dir.create( outputDir , showWarnings = TRUE, recursive = TRUE) 
    }
    
    allGeneSet <- c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA")
    
    noGeneSet <- setdiff( geneSet , allGeneSet )
    if( length(noGeneSet)>0 & verbose ) 
        message( "Gene sets are not defined: ", noGeneSet )
    
    for( gs in intersect(geneSet,allGeneSet) )
    {
        if( verbose ) message('\n\n*** Checking ',gs,'...')
        data(list=gs)
        net <- get(gs)$net
        net <- net[ rownames(net) %in% as.character(data[,geneCol]) , ]
        imputeFC <- data[ match(rownames(net),as.character(data[,geneCol])), ]
        fc <- imputeFC[,fcCol]
        
        if( is.null(weightCol) )
        {
            weights <- rep(1, nrow(net))
        } else {
            weights <- imputeFC[,weightCol]
        }
        
        if( verbose ) message('--> performing SGSEA ...')
        SGSEA.res <- permutationSimpleLmMatrix( fc , net , weights , 
                                              permutationNum, verbose=verbose)
        if( !is.null(get(gs)$annot) )
        {
            annot <- get(gs)$annot
            if( all(table(annot[,1])==1) )
            SGSEA.res <- merge( annot , SGSEA.res , by.x=colnames(annot)[1] ,
                                by.y=colnames(SGSEA.res)[1] )
        }
        SGSEA.res <- SGSEA.res[order(SGSEA.res$empiricalPval) , ]
        write.table( SGSEA.res , paste0(outputDir,'/',gs,'.SGSEA.txt') , 
                     sep='\t' , quote=FALSE , row.names=FALSE , col.names=TRUE)
        
        if( !is.null(MGSEAthres) )
        {
            if( ncol(net)<MGSEAthres )
            {
                if( verbose ) message('\n--> performing MGSEA ...')
                MGSEA.res <- permutationMultipleLmMatrix( fc , net , weights , 
                                              permutationNum, verbose=verbose)
                if( !is.null(get(gs)$annot) )
                {
                    annot <- get(gs)$annot
                    if( all(table(annot[,1])==1) )
                    MGSEA.res <- merge( annot , MGSEA.res , 
                                        by.x=colnames(annot)[1] , 
                                        by.y=colnames(MGSEA.res)[1] )
                }
                MGSEA.res <- MGSEA.res[order(MGSEA.res$empiricalPval) , ]
                write.table( MGSEA.res, paste0(outputDir,'/',gs,'.MGSEA.txt'), 
                             sep='\t' , quote=FALSE , row.names=FALSE , 
                             col.names=TRUE)
            }
        }
    }
    
    return(TRUE)
    
}

