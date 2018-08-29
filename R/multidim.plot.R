#' @export
multidim.plot <- function(files, x.fileNumber=0, y.fileNumber=1) {
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'automation.py', sep='/')
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
