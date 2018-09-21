#' @title GO enrichment of clusters
#' @description perform GO enrichment on cluster groups and plot the result as dot plot
#' @param clusterGroups clustering result from DGE.clust
#' @param OrgDb object returned by AnnotationHub.
#' @param keyType type of genes. e.g. 'ENSEMBL', 'ENTREZID', 'SYMBOL'.
#' @param BgGenes background genes for enrichment analysis
#' @param ont ont One of 'ALL', 'BP', 'BF', or 'CC'.
#' @param top number of top terms to show for each cluster
#' @export
#' @import clusterProfiler
#' @import ggplot2
cluster.enrich <- function(clusterGroups, OrgDb, keyType, BgGenes, ont='BP', top=10){
  # GO enrichment
  universe <- sapply(BgGenes, as.character)
  GO.cluster.res <- as.data.frame(compareCluster(clusterGroups, fun='enrichGO', OrgDb=org.db, keyType=keyType, ont=ont, universe=universe))

  # select top significant terms
  gene.numbers <- c()
  for (n in 1:length(clusterGroups)) {
    gene.numbers <- c(gene.numbers, paste('(', length(clusterGroups[[n]]), ')', sep=''))
  }
  temp <- split(GO.cluster.res, GO.cluster.res$Cluster)
  temp2 <- data.frame(matrix(ncol=ncol(GO.cluster.res), nrow=0))
  colnames(temp2) <- colnames(GO.cluster.res)
  for (n in 1:length(temp)){
    # reorder each sublist by p.adjust and select first top terms
    cluster <- temp[[n]]
    if (nrow(cluster) == 0)
      next
    cluster <- cluster[order(cluster[, 'p.adjust']),]
    if (nrow(cluster) > top){
      cluster <- cluster[1:top,]
    }
    cluster['Cluster'] <- lapply(cluster['Cluster'], function(x) paste(x, gene.numbers[n]))
    temp2 <- rbind(temp2, cluster)
  }
  GO.cluster.res <- temp2[order(temp2[, 'Cluster']),]

  # dot plotting
  GO.cluster.res['Percentage'] <- sapply(GO.cluster.res$GeneRatio, function(x) round(eval(parse(text=x)), 2))
  GO.cluster.res['Description'] <- factor(GO.cluster.res$Description, rev(as.character(unique(GO.cluster.res$Description))))
  ggplot(GO.cluster.res, aes(x=Cluster, y=Description, size=Percentage)) +
    geom_point() +
    theme_bw() +
    aes_string(color='p.adjust') +
    scale_color_continuous(low="red", high='blue', guide=guide_colorbar(reverse=TRUE)) +
    theme(plot.title=element_text(size=14, face='bold', hjust=0.5),
          legend.title=element_text(size=16),
          axis.title=element_text(size=16),
          axis.title.x=element_text(vjust=-0.5, face='bold'),
          axis.title.y=element_text(vjust=2.5, face='bold'),
          axis.text.x=element_text(angle=45, hjust=1, size=10, face='bold'),
          axis.text.y=element_text(size=10, face='bold'))
}

