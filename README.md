# GIGSEA
Genotype Imputed Gene Set Enrichment Analysis

## Description
Various methods of gene set analysis for trait-associated SNPs have been proposed, however, many challenges and limitations remained: 
1. Gene boundaries: different criteria have been proposed to assign a SNP to a gene but no consensus was reached; 
2. Long-range regulation: assigning a causal link to the gene nearest the associated variant falls short of elucidating a long-range functional connection; 
3. Gene size: longer genes are more likely to have significant P-values, possibly inflating the association test for gene sets that have many long genes; 
4. Multiple-marker regulation: the best strategy is not determined on the number of SNPs for each gene and aggregation of different effect sizes of SNPs; 
5. Linkage disequilibrium (LD): the local LD may reduce power to detect associations dependent on multiple markers. 
6. Redundancy among gene sets: a gene may function in multiple ways and thus appear multiple times in functional gene sets. In spite of reflecting the crosstalk between gene sets, the overlap in gene sets may make the results of gene set enrichment analysis more difficult to interpret; 
7. Permutation efficiency: the computational burden of permutation can be substantial; 
8. Threshold-selection: a threshold-dependent procedure may cause the instability of results. 

Here, we presented the GSEA-GS (Gene Set Enrichment Analysis on GWAS Summary Level Data), a novel method that uses GWAS summary statistics and eQTL to infer differential gene expression and interrogate gene set enrichment for the trait-associated SNPs. By incorporating empirical eQTL of disease-relevant tissue, GSEA-GS naturally accounts for factors such as gene size, gene boundary, SNP distal regulation, and multiple-marker regulation. The weighted linear regression model was used to perform the enrichment test, properly adjusting imputation accuracy, model incompleteness and redundancy in different gene sets. The significance level of enrichment is assessed by permutation, where matrix operation was employed to dramatically improve time efficiency. We have shown GSEA-GS has appropriate type I error, and demonstrated high computational efficiency on real data set and discovered the plausible biological findings. 


## Dependencies
-  Matrix

## Installation:
1. Install the [MetaXcan](https://github.com/hakyimlab/MetaXcan) package in Python
2. Install the [devtools](https://github.com/hadley/devtools) package in R
```
   install.packages("devtools")
```
3. Load the devtools package in R
```
   library(devtools) 
```
4. Install GIGSEA in R
```
   install_github("zhushijia/GIGSEA")
```

## Example:
  See [GIGSEA_tutorial](https://github.com/zhushijia/GIGSEA/blob/master/vignettes/GIGSEA_tutorial.Rmd)
