#' runGIGSEA
#'
#' runGIGSEA use MetaXcan to impute the trait-associated differential gene 
#' expression from GWAS summary and eQTL database first, and next, performs 
#' gene set enrichment analysis for the trait-associated SNPs. 
#'
#' @param MetaXcan a character value indicating the path to the MetaXcan.py 
#' file.
#' @param model_db_path a character value indicating the path to tissue 
#' transriptome model.
#' @param covariance a character value indicating the path to file containing 
#' covariance information. This covariance should have information related to 
#' the tissue transcriptome model.
#' @param gwas_folder a character value indicating the folder containing GWAS 
#' summary statistics data.
#' @param gwas_file_pattern a regular expression indicating the gwas summary 
#' files.
#' @param snp_column a character value indicating the name of column holding 
#' SNP data, by default, "SNP".
#' @param non_effect_allele_column a character value indicating the name of 
#' column holding "other/non effect" allele data, by default, "A2".
#' @param effect_allele_column a character value indicating the name of column 
#' holding effect allele data, by default, "A1".
#' @param or_column a character value indicating the name of column holding Odd 
#' Ratio data, by default, "OR".
#' @param beta_column a character value indicating the name of column holding 
#' beta data, by default, "BETA".
#' @param beta_sign_column a character value indicating the name of column 
#' holding sign of beta, by default, "direction".
#' @param zscore_column a character value indicating the name of column holding 
#' zscore of beta, by default, "Z".
#' @param pvalue_column a character value indicating the name of column holding 
#' p-values data, by default, "P".
#' @param gene_set a vector of characters indicating the gene sets of interest 
#' for enrichment test, by default, c("MSigDB.KEGG.Pathway","MSigDB.TF",
#' "MSigDB.miRNA","Fantom5.TF","TargetScan.miRNA","GO", "LINCS.CMap.drug")
#' @param permutation_num an integer indicating the number of permutation.
#' @param output_dir a character value indicating the directory for saving the 
#' results.
#' @param MGSEA_thres an integer value indicating the thresfold for performing 
#' MGSEA. When the number of gene sets is smaller than MGSEAthres, we perform 
#' MGSEA.
#' @param verbose an boolean value indicating whether or not to print output to 
#' the screen 
#'
#' @return TRUE
#' @export
#'
#' @examples
#' # runGIGSEA( MetaXcan="/MetaXcan/software/MetaXcan.py" , 
#'   # model_db_path="data/DGN-WB_0.5.db" ,
#'   # covariance="data/covariance.DGN-WB_0.5.txt.gz" ,
#'   # gwas_folder="data/GWAS" ,
#'   # gwas_file_pattern="heart.summary" ,
#'   # zscore_column="Z" ,
#'   # output_dir="./GIGSEA",
#'   # permutation_num=1000)
#'  
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references 
#' Barbeira, A., et al. Integrating tissue specific mechanisms into GWAS 
#' summary results. bioRxiv 2016:045260.
#' https://github.com/hakyimlab/MetaXcan
#' 
#' @seealso \code{\link{weightedGSEA}}; 
#' 
runGIGSEA <- function( MetaXcan, 
                       model_db_path, 
                       covariance, 
                       gwas_folder, 
                       gwas_file_pattern, 
                       snp_column="SNP", 
                       non_effect_allele_column="A2", 
                       effect_allele_column="A1", 
                       or_column="OR", 
                       beta_column="BETA", 
                       beta_sign_column="direction", 
                       zscore_column="Z", 
                       pvalue_column="P", 
                       gene_set=c("MSigDB.KEGG.Pathway","MSigDB.TF",
                                  "MSigDB.miRNA","TargetScan.miRNA"), 
                       permutation_num=1000, 
                       output_dir="./GIGSEA", 
                       MGSEA_thres=NULL ,
                       verbose = TRUE )
{
    
    MetaXcanCmd<-paste0( MetaXcan , 
                       " --model_db_path " , model_db_path , 
                       " --covariance " , covariance ,
                       " --gwas_folder " , gwas_folder ,
                       " --gwas_file_pattern " , gwas_file_pattern ,
                       " --snp_column " , snp_column ,
                       " --non_effect_allele_column ", non_effect_allele_column,
                       " --effect_allele_column " , effect_allele_column ,
                       " --or_column " , or_column ,
                       " --beta_column " , beta_column ,
                       " --beta_sign_column " , beta_sign_column ,
                       " --zscore_column " , zscore_column ,
                       " --pvalue_column " , pvalue_column ,
                       " --output_file " , output_dir , "/MetaXcan.res.csv" ) 
      
    if( !file.exists(output_dir) )
    { 
        if(verbose) message('creating ' , output_dir )
        dir.create( output_dir , showWarnings = TRUE, recursive = TRUE) 
    }
    
    if(verbose) message(MetaXcanCmd)
    system2(MetaXcanCmd)
    
    metaXcan <- read.table( paste0(output_dir,"/MetaXcan.res.csv") , 
                              sep=',' , header=TRUE )
    cc <- table( as.character( metaXcan$gene_name ) )
    metaXcan <- subset( metaXcan , 
               ! as.character(metaXcan$gene_name) %in% (names(cc)[cc>1]) )
      
    gene <- metaXcan$gene_name
    fc <- metaXcan$zscore
    usedFrac <- metaXcan$n_snps_used / metaXcan$n_snps_in_cov
    r2 <- metaXcan$pred_perf_r2
    weights <- usedFrac*r2
    data <- data.frame(gene,fc,weights)
    
    if(verbose) message( 'permutation for' , permutation_num , 'times' )
    if(verbose) message( 'writing results at' , output_dir )
    
    weightedGSEA(data, geneCol='gene', fcCol='fc', weightCol= 'weights',
                 geneSet=gene_set, permutationNum=permutation_num, 
                 outputDir=output_dir, MGSEAthres=MGSEA_thres)
    
    return(TRUE)
}

