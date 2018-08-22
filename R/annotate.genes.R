#' @export
#' @import clusterProfiler
annotate.genes <- function(OrgDb, keyType, pairedDataset, geneCol, GOEnrichment=FALSE, BgGenes=NULL, ont='ALL'){
  dataset <- pairedDataset
  universe <- BgGenes
  ds.keys <- sapply(dataset[paste(geneCol, 'x', sep='_')], as.character)
  unv.keys <- sapply(universe, as.character)

  if (GOEnrichment == FALSE){
    # assign GO terms directly
    GO.res <- select(OrgDb, ds.keys, c(keyType, 'GO'), keyType)
    GO.res <- GO.res[complete.cases(GO.res['GO']),]
  }
  else {
    # GO enrichment
    GO.res <- enrichGO(gene=dataset$EntrezID, OrgDb=OrgDb, keyType=keyType, ont=ont, universe=universe)
  }

  # construct binary GO term matrix
  for (row in 1:nrow(GO.res)){
    GO.term <- GO.res[row, 'GO']
    if (!GO.term %in% colnames(dataset))
      dataset[GO.term] <- 0
    dataset[dataset[paste(colname, 'x', sep='_')]==GO.res[row, 'SYMBOL'], GO.term] <- 1
  }
  return(dataset)
}
