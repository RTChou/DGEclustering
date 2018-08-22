#' @export
automation <- function(rootDir, geneCol, qvalue=0.05){
  path <- paste(system.file(package="DGEclustering"), "automation.py", sep="/")
  system(paste(path,
               '-d', roodDir,
               '-g', geneCol,
               '-q', qvalue))
}
