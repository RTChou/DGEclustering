#' @export
automation <- function(rootDir, geneCol, qvalue_x=0.05, qvalue_y=0.05, qq_plot=TRUE, fish_plot=TRUE, scatter_plot=TRUE){
  path <- paste(system.file(package="DGEclustering"), "automation.py", sep="/")
  system(paste(path,
               '-d', roodDir,
               '-g', geneCol,
               '-x', qvalue_x,
               '-y', qvalue_y,
               '-q', qq_plot,
               '-f', fish_plot,
               '-s', scatter_plot))
}
