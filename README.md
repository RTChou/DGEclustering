DGEclustering
<img src="../assets/logo.png" height="180" align="right" />
=============
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
Alternatively, the user can install DGEclustering from bioconda, but it may not be the latest version:<br>
conda install r-dgeclustering

## Usage
1. Automation of diagnostic plots: `automation`
   This function will create plotting folders as well as a folder containing paired DGE datasets (merged): `qq_plots`, 
   `fish_plots`, `scatter_plots`, and `paired_files` 
``` r
# Specify where to store the result
dir <- './'

# Automated plotting of diagnostic plots
automation(rootDir=dir, geneCol=gene.col, x.threshold=x.threshold, y.threshold=y.threshold, adjPvalue=adjPvalue, qqPlot=TRUE, fishPlot=TRUE, scatterPlot=TRUE)

# (optional) search for files in the database
mydb <- dbConnect(RSQLite::SQLite(), 'rnaseq.db')
dbListTables(mydb)

# remember to close and disconnect the database
dbDisconnect(mydb)
unlink('rnaseq.db')
```
