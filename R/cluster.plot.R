#' @export
cluster.plot <- function(outDir, filePath1, filePath2, clusterFilePath, geneCol, sigData='ALL', qvalue=0.05, color='brg'){
  path <- paste(system.file(package="DGEclustering"), "cluster_plot.py", sep="/")
  system(paste(path,
               '-d', outDir,
               '-f1', filePath1,
               '-f2', filePath2,
               '-r', clusterFilePath,
               '-g', geneCol,
               '-s', sigData,
               '-q', qvalue,
               '-c', color))
}
