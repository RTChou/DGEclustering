DGEclustering
=============
<img src="../assets/logo.png" height="180" align="right" />
DGEclustering is an R package for multidimensional clustering of differential gene expression (DGE) datasets, and it integrates GO annotations to improve the clustering result.

## Introduction
DGEclustering performs two primary tasks:
1. Searching from a starting directory for differential gene expression datasets (DESeq2 outputs), and automatically generating 
   diagnostic plots for paired DGE datasets, which includes: Q-Q plots, fish plots, and scatter plots.
   ![](../assets/automation.png)<!-- -->
   
2. Conducting multidimensional clustering analysis integrated with GO terms for paired DGE datasets. 

## Installation
1. Create conda environment:<br>
conda create -n r-dgeclustering --file DGEclustering/requirement.txt

2. Install package in R:
``` r
install.packages('DGEclustering', type='source', repos=NULL)
```

