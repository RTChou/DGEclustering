#' @export
cluster.bar.plot <- function(method1.groups, method2.groups, method1.name='method 1', method2.name='method 2', color=NULL){
  percent.matrix <- matrix(NA, nrow=length(method1.groups), ncol=length(method2.groups))
  for (i in 1:length(method1.groups)){
    for (j in 1:length(method2.groups))
      percent.matrix[i, j] <- sum(unlist(method2.groups[j]) %in% unlist(method1.groups[i]))
  }
  stackedMatrix <- data.frame(percent.matrix)
  colnames(stackedMatrix) <- paste('cluster', 1:ncol(percent.matrix))
  stackedMatrix <- stack(stackedMatrix)
  stackedMatrix[, 'ind2'] <- rep(paste('cluster', 1:nrow(percent.matrix)), ncol(percent.matrix))

  if (is.null(color))
    color <- c('#0000ff', '#4000bf', '#80007f', '#c0003f', '#fe0100', '#be4100', '#7e8100', '#3ec100', '#00ff00')

  ggplot(data=stackedMatrix, aes(x=ind2, y=values, fill=ind)) +
    geom_bar(stat="identity", alpha=0.6) +
    theme_classic() +
    scale_fill_manual(values=color) +
    labs(title=(paste('Comparison of Clustering Algorithms\n(', method1.name, ' vs ', method2.name, ')', sep='')))+
    xlab(paste(method1.name, 'clusters')) +
    ylab('number of genes') +
    guides(fill=guide_legend(title=paste(method2.name, 'clusters', sep=' '))) +
    theme(plot.title=element_text(size=14, face='bold', hjust=0.5),
          legend.title=element_text(size=12),
          axis.title=element_text(size=12),
          axis.text.x=element_text(angle=45, hjust=1, color='black'),
          axis.text.y=element_text(color='black'))
}
