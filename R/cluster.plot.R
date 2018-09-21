#' @title visualizating the clustering result
#' @description given the list of dataset file paths or the file path to the MCA result, 
#'   and also the path to the clustering result file, 
#'   this will generate a set clustered scatter/MCA plots.
#' @param outDir directory for storing the output plots
#' @param filePaths a list of dataset file paths or the file path to the MCA result
#' @param MCA If FALSE (which is Default), generate clustered scatter plots
#' @param x.fileNumber the index of the file in the list to be plotted as x axis
#' @param y.fileNumber the index of the file in the list to be plotted as y axis
#' @param geneCol name of the column where the genes are stored in the datasets
#' @param clusterFilePath the file path to the clustering result
#' @param adjPvalue if TRUE, use adjusted p-value for the threshold; if FALSE, use p-value
#' @param color type of matplotlib colormap
#' @export
#' @examples  \dontrun{}
cluster.plot <- function(outDir='./', filePaths, MCA=FALSE, x.fileNumber=1, y.fileNumber=2, geneCol, clusterFilePath, adjPvalue=TRUE, color='brg'){
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'cluster_plot.py', sep='/')
  system(paste(path,
               '-d', outDir,
               '-f', paste(filePaths, collapse=' '),
               '-m', python.boolean.convert(MCA),
               '-n1', x.fileNumber - 1,
	       '-n2', y.fileNumber - 1,
	       '-g', geneCol,
	       '-r', clusterFilePath,
               '-a', python.boolean.convert(adjPvalue),
               '-c', color))
}

