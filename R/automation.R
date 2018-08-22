#' @export
automation <- function(rootDir, geneCol, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE, qqPlot=TRUE, fishPlot=TRUE, scatterPlot=TRUE){
  path <- paste(system.file(package="DGEclustering"), "automation.py", sep="/")
  system(paste(path,
               '-d', roodDir,
               '-g', geneCol,
               '-x', x.threshold,
               '-y', y.threshold,
               '-a', adjPvalue,
               '-q', qqPlot,
               '-f', fishPlot,
               '-s', scatterPlot))
}
