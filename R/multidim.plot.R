#' @export
multidim.plot <- function(filePaths, x.fileNumber=0, y.fileNumber=1, plotDir, datDir, outputName, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE) {
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'automation.py', sep='/')
  system(paste(path,
               '-f', paste(filePaths, collapse=' '),
               '-n1', x.fileNumber,
               '-n2', y.fileNumber,
               '-p', plotDir,
               '-d', datDir,
               '-o', outputName,
               '-x', x.threshold,
               '-y', y.threshold,
	       '-y', python.boolean.convert(adjPvalue)))
}
