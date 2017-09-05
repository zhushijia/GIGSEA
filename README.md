# GIGSEA
Genotype Imputed Gene Set Enrichment Analysis

## Description
We presented the Genotype-imputed Gene Set Enrichment Analysis (GIGSEA), a novel method that uses GWAS-and-eQTL-imputed trait-associated differential gene expression to interrogate gene set enrichment for the trait-associated SNPs. By incorporating eQTL from large gene expression studies, e.g. GTEx, GIGSEA appropriately addresses such challenges for SNP enrichment as gene size, gene boundary, SNP distal regulation, and multiple-marker regulation. The weighted linear regression model, taking as weights both imputation accuracy and model completeness, was used to perform the enrichment test, properly adjusting the bias due to redundancy in different gene sets. The permutation test, furthermore, is used to evaluate the significance of enrichment, whose efficiency can be largely elevated by expressing the computational intensive part in terms of large matrix operation. We have shown the appropriate type I error rates for GIGSEA (<5%), and the preliminary results also demonstrate its good performance to uncover the real signal. 

## Dependencies
-  Matrix

## Installation:
1. Install the [MetaXcan](https://github.com/hakyimlab/MetaXcan) package
2. Install the [devtools](https://github.com/hadley/devtools) package
```
   install.packages("devtools")
```
3. Load the devtools package
```
   library(devtools)
```
4. Install GIGSEA
```
   install_github("zhushijia/GIGSEA")
```

## Example:
  See [GIGSEA_tutorial](https://github.com/zhushijia/GIGSEA/blob/master/vignettes/GIGSEA_tutorial.Rmd)
