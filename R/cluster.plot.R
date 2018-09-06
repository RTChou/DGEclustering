#' @export
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

