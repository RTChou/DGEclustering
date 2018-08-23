#' @export
automation <- function(rootDir, geneCol, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE, qqPlot=TRUE, fishPlot=TRUE, scatterPlot=TRUE){
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package="DGEclustering"), "automation.py", sep="/")
  system(paste(path,
               '-d', rootDir,
               '-g', geneCol,
               '-x', x.threshold,
               '-y', y.threshold,
               '-a', python.boolean.convert(adjPvalue),
               '-q', python.boolean.convert(qqPlot),
               '-f', python.boolean.convert(fishPlot),
               '-s', python.boolean.convert(scatterPlot)))
}
