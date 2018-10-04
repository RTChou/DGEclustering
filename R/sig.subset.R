#' @export
#' @import png
sig.subset <- function(datasets, geneCol, x.fileNumber=1, y.fileNumber=2, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE) {
  # export datasets to temp folder
  tempFolder <- '/tmp/dgeclustering'
  system(paste('mkdir -p', tempFolder))
  file.path(tempFolder, paste0('ds', seq(1,length(datasets)), '.tsv'))
  for (i in 1:length(datasets)){ 
    write.table(datasets[i], file=filepaths[i], sep='\t', row.names=FALSE)
  }
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'sig_subset.py', sep='/')
  system(paste(path,
             '-f', paste(filePaths, collapse=' '),
             '-g', geneCol,
	     '-n1', x.fileNumber - 1,
             '-n2', y.fileNumber - 1,
             '-p', tempFolder,
             '-d', tempFolder,
             '-o', 'temp',
             '-x', x.threshold,
             '-y', y.threshold,
             '-y', python.boolean.convert(adjPvalue)))
}

