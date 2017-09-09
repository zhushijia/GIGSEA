<<<<<<< HEAD
#' runGIGSEA
#'
#' runGIGSEA use MetaXcan to impute the trait-associated differential gene expression from GWAS summary and eQTL database first, and next, performs gene set enrichment analysis for the trait-associated SNPs. 
#'
#' @param MetaXcan a character value indicating the path to the MetaXcan.py file.
#' @param model_db_path a character value indicating the path to tissue transriptome model.
#' @param covariance a character value indicating the path to file containing covariance information. This covariance should have information related to the tissue transcriptome model.
#' @param gwas_folder a character value indicating the folder containing GWAS summary statistics data.
#' @param gwas_file_pattern a regular expression indicating the gwas summary files.
#' @param snp_column a character value indicating the name of column holding SNP data, by default, "SNP".
#' @param non_effect_allele_column a character value indicating the name of column holding "other/non effect" allele data, by default, "A2".
#' @param effect_allele_column a character value indicating the name of column holding effect allele data, by default, "A1".
#' @param or_column a character value indicating the name of column holding Odd Ratio data, by default, "OR".
#' @param beta_column a character value indicating the name of column holding beta data, by default, "BETA".
#' @param beta_sign_column a character value indicating the name of column holding sign of beta, by default, "direction".
#' @param zscore_column a character value indicating the name of column holding zscore of beta, by default, "Z".
#' @param pvalue_column a character value indicating the name of column holding p-values data, by default, "P".
#' @param gene_set a vector of characters indicating the gene sets of interest for enrichment test, by default, c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA","Fantom5.TF","TargetScan.miRNA","GO", "LINCS.CMap.drug")
#' @param permutation_num an integer indicating the number of permutation.
#' @param output_dir a character value indicating the directory for saving the results.
#' @param MGSEA_thres an integer value indicating the thresfold for performing MGSEA. When the number of gene sets is smaller than MGSEAthres, we perform MGSEA.
#'
#' @return
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
#' Barbeira, A., et al. Integrating tissue specific mechanisms into GWAS summary results. bioRxiv 2016:045260.
#' https://github.com/hakyimlab/MetaXcan
#' 
#' @seealso \code{\link{GIGSEA}}; 
#' 
runGIGSEA <- function( MetaXcan , model_db_path, covariance, gwas_folder, gwas_file_pattern, 
                       snp_column="SNP", non_effect_allele_column="A2", effect_allele_column="A1", or_column="OR", beta_column="BETA", beta_sign_column="direction", zscore_column="Z", pvalue_column="P", 
                       gene_set=c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA","Fantom5.TF","TargetScan.miRNA","GO","LINCS.CMap.drug"), permutation_num=1000, output_dir="./GIGSEA", MGSEA_thres=NULL )
{

  MetaXcanCmd = paste0( MetaXcan , 
	  " --model_db_path " , model_db_path , 
		" --covariance " , covariance ,
		" --gwas_folder " , gwas_folder ,
		" --gwas_file_pattern " , gwas_file_pattern ,
		" --snp_column " , snp_column ,
		" --non_effect_allele_column " , non_effect_allele_column ,
		" --effect_allele_column " , effect_allele_column ,
		" --or_column " , or_column ,
		" --beta_column " , beta_column ,
		" --beta_sign_column " , beta_sign_column ,
		" --zscore_column " , zscore_column ,
		" --pvalue_column " , pvalue_column ,
		" --output_file " , output_dir , "/MetaXcan.res.csv" ) 
	
	dir.create( output_dir , showWarnings = TRUE, recursive = TRUE)
	preWD = setwd( output_dir )
	
	cat(MetaXcanCmd,'\n')
	system(MetaXcanCmd)

	
	metaXcan <- read.table( paste0(output_dir,"/MetaXcan.res.csv") , sep=',' , header=T )
	gene <- metaXcan$gene_name
	fc <- metaXcan$zscore
	usedFrac <- metaXcan$n_snps_used / metaXcan$n_snps_in_cov
	r2 <- metaXcan$pred_perf_r2
	weights <- usedFrac*r2
	data <- data.frame(gene,fc,weights)

	cat( 'permutation for' , permutation_num , 'times\n' )
	cat( 'writing results at' , output_dir , '\n' )

	GIGSEA(data, geneCol='gene', fcCol='fc', weightCol= 'weights', 	geneSet=gene_set, permutationNum=permutation_num, outputDir=output_dir, MGSEAthres=MGSEA_thres)

	setwd(preWD)
	
}

# source ~/setup/python/virtualenv-1.10.1/myVE/bin/activate
# module load py_packages

# runGIGSEA( MetaXcan="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/PrediXcan/MetaXcan/software/MetaXcan.py" , 
	# model_db_path="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/PrediXcan/MetaXcan/software/data/DGN-WB_0.5.db" ,
	# covariance="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/PrediXcan/MetaXcan/software/data/covariance.DGN-WB_0.5.txt.gz" ,
	# gwas_folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/TWAS/source/fusion_twas-master" ,
	# gwas_file_pattern="PGC2.SCZ.sumstats" ,
	# zscore_column="Z" ,
	# output_dir="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/GIGSEA",
	# permutation_num=1000)

# library(GIGSEA , lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/tmp")

# MetaXcan="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/PrediXcan/MetaXcan/software/MetaXcan.py" 
# model_db_path="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/PrediXcan/MetaXcan/software/data/DGN-WB_0.5.db" 
# covariance="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/PrediXcan/MetaXcan/software/data/covariance.DGN-WB_0.5.txt.gz" 
# gwas_folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/TWAS/source/fusion_twas-master" 
# gwas_file_pattern="PGC2.SCZ.sumstats"
# zscore_column="Z" 
# gene_set=c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA","Fantom5.TF","TargetScan.miRNA","GO","LINCS.CMap.drug")
# geneSet=c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA","Fantom5.TF","TargetScan.miRNA","GO","LINCS.CMap.drug")
# output_dir="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/GIGSEA"
# permutation_num=1000

