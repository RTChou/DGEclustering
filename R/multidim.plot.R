#' @export
multidim.plot <- function(filePaths, x.fileNumber=1, y.fileNumber=2, plotDir, datDir, outputName=NULL, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE) {
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'automation.py', sep='/')
  if (!is.null(outputName)) {
    system(paste(path,
               '-f', paste(filePaths, collapse=' '),
               '-n1', x.fileNumber - 1,
               '-n2', y.fileNumber - 1,
               '-p', plotDir,
               '-d', datDir,
               '-o', outputName,
               '-x', x.threshold,
               '-y', y.threshold,
               '-y', python.boolean.convert(adjPvalue)))
  }
  else{
    system(paste(path,
               '-f', paste(filePaths, collapse=' '),
               '-n1', x.fileNumber - 1,
               '-n2', y.fileNumber - 1,
               '-p', plotDir,
               '-d', datDir,
               '-x', x.threshold,
               '-y', y.threshold,
               '-y', python.boolean.convert(adjPvalue)))
  }
}

