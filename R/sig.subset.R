#' @title
#' @description
#' @param datasets
#' @param geneCol the column name of gene IDs
#' @param x.fileNumber
#' @param y.fileNumber
#' @param x.threshold (adj)p-value threshold of the first dataset
#' @param y.threshold (adj)p-value threshold of the second dataset
#' @param adjPvalue if TRUE, use adjusted p-value for the threshold; if FALSE, use p-value
#' @export
#' @import png
sig.subset <- function(datasets, geneCol, x.fileNumber=1, y.fileNumber=2, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE) {
  # export datasets to temp folder
  temp.folder <- '/tmp/dgeclustering'
  system(paste('mkdir -p', temp.folder))
  filepaths <- file.path(temp.folder, names(datasets), '.tsv')
  for (i in 1:length(datasets)){ 
    write.table(datasets[i], file=filepaths[i], sep='\t', row.names=FALSE)
  }
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'sig_subset.py', sep='/')
  system(paste(path,
             '-f', paste(filepaths, collapse=' '),
             '-g', geneCol,
	     '-n1', x.fileNumber - 1,
             '-n2', y.fileNumber - 1,
             '-p', temp.folder,
             '-d', temp.folder,
             '-o', 'temp',
             '-x', x.threshold,
             '-y', y.threshold,
             '-y', python.boolean.convert(adjPvalue)))
  # plot the scatter plot
  img <- readPNG(file.path(temp.folder, 'temp_sig_plot.png'))
  plot.new() 
  rasterImage(img,0,0,1,1)
  if (length(datasets) == 2) {
    dis <- read.table(file.path(temp.folder, 'temp_disagreeing_genes.tsv'), header=TRUE, 
      check.names=FALSE, sep='\t', stringsAsFactors=FALSE)
    con <- read.table(file.path(temp.folder, 'temp_agreeing_genes.tsv'), header=TRUE, 
      check.names=FALSE, sep='\t', stringsAsFactors=FALSE)
    dat <- list(dis, con)
    names(dat) <- c('dis', 'con')
    return(dat)
  } 
  else {
    dat <- read.table(file.path(temp.folder, 'temp_all_sig_genes.tsv'), header=TRUE, check.names=FALSE, 
      sep='\t', stringsAsFactors=FALSE)
    return(dat)
  }
}
