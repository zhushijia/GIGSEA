runGIGSEA <- function( MetaXcan , model_db_path, covariance, gwas_folder, gwas_file_pattern, zscore_column, output_dir, 
	geneSet=c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA","Fantom5.TF","TargetScan.miRNA","GO","LINCS.CMap.drug"), permutation_num=1000 )
{

	cmd = paste0( MetaXcan , " --model_db_path " , model_db_path , 
		" --covariance " , covariance ,
		" --gwas_folder " , gwas_folder ,
		" --gwas_file_pattern " , gwas_file_pattern ,
		" --zscore_column " , zscore_column ,
		" --output_file " , output_dir , "/MetaXcan.res.csv" ) 
	
	dir.create( output_dir , showWarnings = TRUE, recursive = TRUE)
	preWD = setwd( output_dir )
	
	cat(cmd,'\n')
	system(cmd)
	
	metaXcan <- read.table( paste0(output_dir,"/MetaXcan.res.csv") , sep=',' , header=T )
	gene <- metaXcan$gene_name
	fc <- metaXcan$zscore
	usedFrac <- metaXcan$n_snps_used / metaXcan$n_snps_in_cov
	r2 <- metaXcan$pred_perf_r2
	weights <- usedFrac*r2
	data <- data.frame(gene,fc,weights)

	cat( 'permutation for' , permutation_num , 'times\n' )
	cat( 'writing results at' , output_dir , '\n' )

	GIGSEA(data, geneCol='gene', fcCol='fc', weightCol= 'weights', 	geneSet=geneSet, permutationNum=permutation_num, outputDir=output_dir )

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
# geneSet=c("MSigDB.KEGG.Pathway","MSigDB.TF","MSigDB.miRNA","Fantom5.TF","TargetScan.miRNA","GO","LINCS.CMap.drug")
# output_dir="/hpc/users/zhus02/schzrnas/sjzhu/Project/GWAS/GIGSEA"
# permutation_num=1000

