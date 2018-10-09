DGEclustering
<img src="../assets/logo.png" height="180" align="right" />
=============
DGEclustering is an R package for multidimensional clustering of differential gene expression (DGE) datasets, and it integrates GO annotations to improve the clustering result.


## Introduction
DGEclustering performs two primary tasks:<br>
<b>I.</b> Searching from a starting directory for differential gene expression datasets (DESeq2 outputs), and automatically generating diagnostic plots for paired DGE datasets, which includes: Q-Q plots, fish plots, and scatter plots.
<p align="center"><img src="../assets/automation.png" width="700"></p>
   
<b>II.</b> Conducting multidimensional clustering analysis integrated with GO terms for paired (could be two or more) DGE datasets. The GO terms are intergrated using the InteGO package (Verbanck, M., Lê, S. & Pagès, J. A new unsupervised gene clustering algorithm based on the integration of biological knowledge into expression data. BMC Bioinformatics 14, 42 (2013).)
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
automation(rootDir=dir, geneCol=gene.col, x.threshold=x.threshold, y.threshold=y.threshold, 
adjPvalue=adjPvalue, qqPlot=TRUE, fishPlot=TRUE, scatterPlot=TRUE)

# (optional) search for files in the database
mydb <- dbConnect(RSQLite::SQLite(), 'rnaseq.db')
dbListTables(mydb)

# remember to close and disconnect the database
dbDisconnect(mydb)
unlink('rnaseq.db')
```

### II. Multidimensional clustering integrated with GO terms
#### Step 0-1: Specify the column name for gene IDs (for DGE datasets), organism database, and key type
``` r
gene.col <- 'gene'
orgdb <- query(hub, 'Drosophila melanogaster')[['AH57972']]
keytype <- 'FLYBASE'
```

#### Step 0-2: Specify the paired input files 
``` r
# Import example datasets
data(list=c('treatment1.vs.control', 'treatment2.vs.control', 'treatment3.vs.control'))

# Create a list of datasets
## in this example we cluster on only two paired datasets
datasets <- list(treatment1.vs.control, treatment2.vs.control)
names(datasets) <- c('treatment1.vs.control', 'treatment2.vs.control')
```

#### Step 1-1: Set the (adjusted) p-value thresholds for selecting genes
``` r
x.threshold <- 0.05
y.threshold <- 0.05
adjPvalue <- TRUE
```

#### Step 1-2: Merge datasets and subset the significant genes: `sig.subset`
This function Generates significant scatter plot as well as merged significant subsets
``` r
## 'x.dsNumber' and 'y.dsNumber' is for plotting purposes (x- and y-axis)
sig.res <- sig.subset(datasets, geneCol=gene.col, x.dsNumber=1, y.dsNumber=2, x.threshold=x.threshold, 
y.threshold=y.threshold, adjPvalue=adjPvalue)

## show Plot
sig.res$p

## show datasets
## for two paired datasets, there will be discordant and concordant datasets. Discordant dataset contains 
## a set of genes having different signs of log2 fold changes between the paired datasets, whereas 
## concordant dataset contains genes having same signs of log2 fold changes.
if (length(sig.res) == 3) { ## for two paired datasets
  head(sig.res$dis)
  head(sig.res$con)
} else { ## For more paired datasets
  head(sig.res$dat)
}
```

#### Step 2: Annotate Genes: `annotate.genes`
``` r
# Assign the `dat` variable
if (length(sig.res) == 3) { ## for two paired datasets we can combine the discordant and condordant datasets
  dat <- rbind(sig.res$dis, sig.res$con)
} else { ## for more paired datasets
  dat <- sig.res$dat
}

# Background genes for GO enrichment
bg.genes <- treatment1.vs.control[gene.col]

# Annotate genes
ann <- annotate.genes(OrgDb=orgdb, keyType=keytype, genes=unlist(dat[gene.col]), 
GOEnrichment=FALSE, BgGenes=bg.genes)
```

#### Step 3: Prepare expression and annotation datasets for clustering
expressions: choose the desired dimensions <br> 
annotations: choose the number of GO terms for optimal clustering <br> 
number of group: choose the number of group for optimal clustering <br> 
``` r
# expression dataset
## show column names for dat
colnames(dat)
## this example selects log2 fold changes and adjusted p-values as input
exp <- dat[,grepl("log2FoldChange|padj", colnames(dat))]
rownames(exp) <- make.names(dat[,gene.col], unique=TRUE)

# annotaiton dataset
## calculate number of GO terms assign to a specific number of genes
sum(apply(ann, 2, sum) >= 100)
## choose the appropriate number of GO terms for annotation dataset
ann <- ann[, apply(ann, 2, sum) >= 100]

# choose the number of clustering groups
nb.group=8
```

#### Step 4: Cluster genes, and visualize the result: `DGE.clust`, `cluster.plot`
``` r
# Clustering analysis
res <- DGE.clust(expressions=exp, annotations=ann, clust.method='intego', nb.group=nb.group)
res$groups ## view clustering result
res$vignette ## view vignette

# visualize the clustering result
p <- cluster.plot(datasets, res$groups, x.dsNumber=1, y.dsNumber=2, geneCol=gene.col, 
adjPvalue=adjPvalue, color='brg')
p ## show plot

# visualize the clustering result (MCA plot)
p.MCA <- cluster.plot(res.groups=res$groups, res.MCA=res$MCA, MCA=TRUE, geneCol=gene.col, 
adjPvalue=adjPvalue, color='brg')
p.MCA ## show plot
```

#### Step 5: GO enrichment of the clustering result (visualization): `cluster.enrich`
``` r
# Background genes for GO enrichment
bg.genes <- treatment1.vs.control[gene.col]

# Visualizing GO enrichment result, and showing only the top 10 terms for each cluster
p.GO <- cluster.enrich(clusterGroups=res$groups, OrgDb=orgdb, keyType=keytype, BgGenes=bg.genes, 
ont='BP', top=10)

# show plot
p.GO
```

### III. Session information
``` r
sessionInfo()
```

```
## R version 3.4.1 (2017-06-30)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Fedora 28 (Twenty Eight)

## Matrix products: default
## BLAS: /home/guest/miniconda3/envs/r-dgeclustering/lib/R/lib/libRblas.so
## LAPACK: /home/guest/miniconda3/envs/r-dgeclustering/lib/R/lib/libRlapack.so

## locale:
##  [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C             
##  [3] LC_TIME=en_US.utf8        LC_COLLATE=en_US.utf8    
##  [5] LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8   
##  [7] LC_PAPER=en_US.utf8       LC_NAME=C                
##  [9] LC_ADDRESS=C              LC_TELEPHONE=C           
## [11] LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C      

## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     

## other attached packages:
##  [1] InteGO_2.0            AnnotationDbi_1.40.0  IRanges_2.12.0       
##  [4] S4Vectors_0.16.0      Biobase_2.38.0        clusterProfiler_3.6.0
##  [7] DOSE_3.4.0            ggplot2_2.2.1         stringr_1.3.0        
## [10] RSQLite_2.0           AnnotationHub_2.10.1  BiocGenerics_0.24.0  
## [13] DGEclustering_0.1.0   cluster_2.0.6         FactoMineR_1.41      

## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.15                  lattice_0.20-34              
##  [3] tidyr_0.8.1                   GO.db_3.5.0                  
##  [5] png_0.1-7                     digest_0.6.15                
##  [7] mime_0.5                      R6_2.2.2                     
##  [9] plyr_1.8.4                    BiocInstaller_1.28.0         
## [11] httr_1.3.1                    pillar_1.2.2                 
## [13] rlang_0.2.1                   curl_3.2                     
## [15] lazyeval_0.2.1                data.table_1.10.4            
## [17] blob_1.1.1                    qvalue_2.10.0                
## [19] splines_3.4.1                 BiocParallel_1.12.0          
## [21] igraph_1.2.1                  bit_1.1-12                   
## [23] munsell_0.5.0                 shiny_1.0.5                  
## [25] fgsea_1.4.0                   compiler_3.4.1               
## [27] httpuv_1.3.6.2                pkgconfig_2.0.1              
## [29] htmltools_0.3.6               flashClust_1.01-2            
## [31] tibble_1.4.2                  gridExtra_2.3                
## [33] interactiveDisplayBase_1.16.0 MASS_7.3-48                  
## [35] leaps_3.0                     grid_3.4.1                   
## [37] xtable_1.8-2                  gtable_0.2.0                 
## [39] DBI_1.0.0                     magrittr_1.5                 
## [41] scales_0.5.0                  rlist_0.4.6.1                
## [43] stringi_1.1.7                 GOSemSim_2.4.0               
## [45] reshape2_1.4.3                scatterplot3d_0.3-41         
## [47] DO.db_2.9                     rvcheck_0.0.9                
## [49] fastmatch_1.1-0               tools_3.4.1                  
## [51] bit64_0.9-5                   glue_1.2.0                   
## [53] purrr_0.2.4                   yaml_2.1.18                  
## [55] colorspace_1.3-2              memoise_1.1.0           
```
