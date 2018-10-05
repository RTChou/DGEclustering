#' @title visualizating the clustering result
#' @description given the list of dataset file paths or the file path to the MCA result, 
#'   and also the path to the clustering result file, 
#'   this will generate a set clustered scatter/MCA plots.
#' @param outDir directory for storing the output plots
#' @param datasets a list of datasets
#' @param MCA If FALSE (which is Default), generate clustered scatter plots
#' @param x.dsNumber the index of the dataset in the list to be plotted as x axis
#' @param y.dsNumber the index of the dataset in the list to be plotted as y axis
#' @param geneCol name of the column where the genes are stored in the datasets
#' @param clusterFilePath the file path to the clustering result
#' @param adjPvalue if TRUE, use adjusted p-value for the threshold; if FALSE, use p-value
#' @param color type of matplotlib colormap
#' @export
#' @examples  \dontrun{}
cluster.plot <- function(datasets, res.groups, res.MCA, MCA=FALSE, x.dsNumber=1, y.dsNumber=2, geneCol, adjPvalue=TRUE, color='brg'){
  temp.folder <- '/tmp/dgeclustering'
  system(paste('mkdir -p', temp.folder))
  if (MCA == FALSE) {
    if (is.null(datasets) | is.null(res.groups)){i
      stop('arguments \"datasets\" and \"res.groups\" should not be NULL.')
    }
    filepaths <- file.path(temp.folder, paste0(names(datasets), '.tsv'))
    cluster.filepath <- file.path(temp.folder, 'clusters.tsv')
    for (i in 1:length(datasets)){
      write.table(datasets[[i]], file=filepaths[i], sep='\t', row.names=FALSE)
    }
  }  
  else {
    if (is.null(res.groups) | is.null(res.MCA)) {
      stop('\"res.groups\" and \"res.MCA\" should not be NULL.')
    }
    filepaths <- file.path(temp.folder, 'MCA.tsv')
    write.table(res$MCA, 'MCA.tsv', sep='\t', row.names=TRUE, col.names=NA)
  }
  write.table(stack(res.groups), file=cluster.filepath, sep='\t', row.names=FALSE)
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'cluster_plot.py', sep='/')
  system(paste(path,
               '-d', temp.folder,
               '-f', paste(filepaths, collapse=' '),
               '-m', python.boolean.convert(MCA),
               '-n1', x.dsNumber - 1,
	       '-n2', y.dsNumber - 1,
	       '-g', geneCol,
	       '-r', cluster.filepath,
               '-a', python.boolean.convert(adjPvalue),
               '-c', color))
  # plotting
  plot <- list()
  img <- readPNG(file.path(temp.folder, 'cluster_all.png'))
  g <- rasterGrob(img, interpolate=TRUE)
  p <- ggplot() +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  list.append(plot, p)
  for (i in (1:length(res.groups))) {
    img <- readPNG(file.path(temp.folder, paste0('cluster_', i, '.png')))
    g <- rasterGrob(img, interpolate=TRUE)
    p <- ggplot() +
      annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
    list.append(plot, p)
  }
  names(plot) <- c('p', paste0('p', seq(1, length(res.groups))))
  return(plot)
}
