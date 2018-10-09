DGEclustering
<img src="../assets/logo.png" height="180" align="right" />
=============
DGEclustering is an R package for multidimensional clustering of differential gene expression (DGE) datasets, and it integrates GO annotations to improve the clustering result.


## Introduction
DGEclustering performs two primary tasks:<br>
### I.  Searching from a starting directory for differential gene expression datasets (DESeq2 outputs), and automatically generating 
    diagnostic plots for paired DGE datasets, which includes: Q-Q plots, fish plots, and scatter plots.
    <p align="center"><img src="../assets/automation.png" width="700"></p>
   
### II. Conducting multidimensional clustering analysis integrated with GO terms for paired DGE datasets. The GO terms are 
    intergrated using the InteGO package (Verbanck, M., Lê, S. & Pagès, J. A new unsupervised gene clustering algorithm based on 
    the integration of biological knowledge into expression data. BMC Bioinformatics 14, 42 (2013).)
    <p align="center"><img src="../assets/clustering.png" width="650"></p>


## Installation
1. Create conda environment:
``` bash
conda create -n r-dgeclustering --file DGEclustering/requirement.txt
```

2. Install package in R:
``` r
install.packages('DGEclustering', type='source', repos=NULL)
```
Alternatively, the user can install DGEclustering from bioconda, but it may not be the latest version:
``` bash
conda install r-dgeclustering
```


## Usage
Import libraries:
``` r
library(DGEclustering)
library(AnnotationHub)
library(RSQLite)
library(stringr)
library(ggplot2)
library(clusterProfiler)
hub <- AnnotationHub::.Hub("AnnotationHub",
        getAnnotationHubOption("URL"),

        # Cache location is specific to this instance of lcdb-wf so we don't
        # clobber old runs with new annotation data
        './AnnotationHubCache',

        httr::use_proxy(Sys.getenv("http_proxy")),
        FALSE)
```

### I. Automation of diagnostic plots: `automation`
   This function will create plotting folders as well as a folder containing paired DGE datasets (merged): `qq_plots`, 
   `fish_plots`, `scatter_plots`, and `paired_files`.
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

### II. Multidimensional clustering integrated with GO terms
#### Step 1: Specify the column name for gene IDs, organism database, and key type
``` r
gene.col <- 'gene'
orgdb <- query(hub, 'Drosophila melanogaster')[['AH57972']]
keytype <- 'FLYBASE'
```

#### Step 2: Set the (adjusted) p-value thresholds
``` r
x.threshold <- 0.05
y.threshold <- 0.05
adjPvalue <- TRUE
```

#### Step 3: Specify the input files and subset significant genes
``` r
# Import example datasets
data(list=c('treatment1.vs.control', 'treatment2.vs.control', 'treatment3.vs.control'))

# Create a list of datasets
datasets <- list(treatment1.vs.control, treatment2.vs.control)
names(datasets) <- c('treatment1.vs.control', 'treatment2.vs.control')

# Generate significant plot and significant subsets
sig.res <- sig.subset(datasets, geneCol=gene.col, x.dsNumber=1, y.dsNumber=2, x.threshold=x.threshold, y.threshold=y.threshold, adjPvalue=adjPvalue)

# show Plot
sig.res$p
```

#### Step 4: Annotate Genes
``` r
## For two paired datasets
if (length(sig.res) == 3) {
  dat <- rbind(sig.res$dis, sig.res$con)
} else { ## For time course data
  dat <- sig.res$dat
}

## Background genes for GO enrichment
bg.genes <- treatment1.vs.control[gene.col]

## Annotate genes
ann <- annotate.genes(OrgDb=orgdb, keyType=keytype, genes=unlist(dat[gene.col]), GOEnrichment=FALSE, BgGenes=bg.genes)
```

#### Step 5: Prepare expression and annotation datasets for clustering
expressions: choose the desired dimensions <br> 
annotations: choose the number of GO terms for optimal clustering <br> 
number of group: choose the number of group for optimal clustering <br> 
``` r
# expressions
exp <- dat[,grepl("log2FoldChange|padj", colnames(dat))]
rownames(exp) <- make.names(dat[,gene.col], unique=TRUE)

# annotaitons
## less than 30 GO terms
sum(apply(ann, 2, sum) >= 5)
ann <- ann[, apply(ann, 2, sum) >= 100]

# number of groups
nb.group=8
```

#### Step 6: Cluster genes, and visualize the result
``` r
# Clustering analysis
res <- DGE.clust(expressions=exp, annotations=ann, clust.method='intego', nb.group=nb.group)

# view clustering result
res$groups

# visualize the clustering result
p <- cluster.plot(datasets, res$groups, x.dsNumber=1, y.dsNumber=2, geneCol=gene.col, adjPvalue=adjPvalue, color='brg')
p

# visualize the clustering result (MCA plot)
p.MCA <- cluster.plot(res.groups=res$groups, res.MCA=res$MCA, MCA=TRUE, geneCol=gene.col, adjPvalue=adjPvalue, color='brg')
p.MCA
```

#### Step 7: GO enrichment of the clustering result (visualization)
``` r
# Background genes for GO enrichment
bg.genes <- read.table(file1, header=TRUE, sep='\t', stringsAsFactors=TRUE)[gene.col]

# visualizing GO enrichment result
p.GO <- cluster.enrich(clusterGroups=res$groups, OrgDb=orgdb, keyType=keytype, BgGenes=bg.genes, ont='BP', top=10)
p.GO
```
