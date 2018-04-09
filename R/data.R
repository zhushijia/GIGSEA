#' heart.metaXcan
#'
#' The MetaXcan-predicted differential gene expression from the cardiovascular disease (CVD) GWAS, CARDIoGRAMplusC4D (60,801 cases, 123,504 controls and 9.4M SNPs). 
#'
#' @format A data frame with the following items:
#' \describe{
#'   \item{gene}{a gene's id}
#'   \item{gene_name}{a gene's name}
#'   \item{zscore}{MetaXcan's association result for the gene}
#'   \item{effect_size}{MetaXcan's association effect size for the gene}
#'   \item{pvalue}{P-value of the aforementioned statistic}
#'   \item{pred_perf_r2}{R2 of transcriptome prediction model's correlation to gene's measured transcriptome}
#'   \item{pred_perf_pval}{pval of transcriptome prediction model's correlation to gene's measured transcriptome}
#'   \item{pred_perf_qval}{qval of transcriptome prediction model's correlation to gene's measured transcriptome}
#'   \item{n_snps_used}{number of snps from GWAS that got used in MetaXcan analysis}
#'   \item{n_snps_in_cov}{number of snps in the covariance matrix}
#'   \item{n_snps_in_model}{number of snps in the prediction model}
#'   \item{var_g}{variance of the gene expression}
#'   ...
#' }
#' @source \url{http://www.cardiogramplusc4d.org/data-downloads/;https://cloud.hakyimlab.org/s-predixcan}
"heart.metaXcan"


#' MSigDB.KEGG.Pathway
#'
#' Gene sets derived from the KEGG pathway database. 
#'
#' @format A list with two items:
#' \describe{
#'   \item{net}{a sparse matrix, the connectivity between terms and genes, comprising 186 pathways (column) and 5267 genes (row)}
#'   \item{annot}{a data frame, description of terms}
#'   ...
#' }
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2}
"MSigDB.KEGG.Pathway"


#' MSigDB.TF
#'
#' Gene sets that share upstream cis-regulatory motifs which can function as potential transcription factor binding sites.
#'
#' @format A list with two items:
#' \describe{
#'   \item{net}{a sparse matrix, the connectivity between terms and genes, comprising 615 TFs (column) and 12774 genes (row)}
#'   \item{annot}{a data frame, description of terms}
#'   ...
#' }
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C3}
"MSigDB.TF"


#' MSigDB.miRNA
#'
#' Gene sets that contain genes sharing putative target sites (seed matches) of human mature miRNA in their 3'-UTRs.
#'
#' @format A list with two items:
#' \describe{
#'   \item{net}{a sparse matrix, the connectivity between terms and genes, comprising 221 miRNAs (column) and 7444 genes (row)}
#'   \item{annot}{a data frame, description of terms}
#'   ...
#' }
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C3}
"MSigDB.miRNA"


#' TargetScan.miRNA
#'
#' Gene sets of predicted human miRNA target sites were downloaded from TargetScan. 
#' TargetScan groups miRNAs that have identical subsequences at positions 2 through 8 of the miRNA, 
#' i.e. the 2-7 seed region plus the 8th nucleotide, and provides predictions for each such seed motif.  
#'
#' @format A list with two items:
#' \describe{
#'   \item{net}{a sparse matrix, the connectivity between terms and genes, comprising 87 miRNA seed motifs and 9861 genes}
#'   \item{annot}{a data frame, description of terms}
#'   ...
#' }
#' @source \url{http://www.targetscan.org}
"TargetScan.miRNA"



