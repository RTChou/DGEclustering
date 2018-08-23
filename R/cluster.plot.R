#' @export
cluster.plot <- function(outDir, filePath1, filePath2, clusterFilePath, geneCol, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE, sigData='ALL', color='brg'){
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package="DGEclustering"), "cluster_plot.py", sep="/")
  system(paste(path,
               '-d', outDir,
               '-f1', filePath1,
               '-f2', filePath2,
               '-r', clusterFilePath,
               '-g', geneCol,
               '-x', x.threshold,
               '-y', y.threshold,
               '-a', python.boolean.convert(adjPvalue),
               '-s', sigData,
               '-c', color))
}
