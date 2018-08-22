#' @export
#' @import AnnotationDbi
#' @import clusterProfiler
#' @import ggplot2
cluster.sim.plot <- function(method1.groups, method2.groups, OrgDb, keytype, ont='BP', method1.name='method 1', method2.name='method 2'){
  sim.matrix <- matrix(NA, nrow=length(method1.groups), ncol=length(method2.groups))
  for (i in 1:length(method1.groups)){
    for (j in 1:length(method2.groups)){
      c1 <- mapIds(OrgDb, keys=unlist(method1.groups[i]), column="ENTREZID", keytype=keytype, multiVals="first")
      c2 <- mapIds(OrgDb, keys=unlist(method2.groups[j]), column="ENTREZID", keytype=keytype, multiVals="first")
      semData <- godata(OrgDb=OrgDb, ont=ont, computeIC=FALSE)
      sim.matrix[i, j] <- clusterSim(c1, c2, semData=semData, measure='Wang')
    }
  }
  rownames(sim.matrix) <- paste('cluster', 1:length(method1.groups))
  colnames(sim.matrix) <- paste('cluster', 1:length(method2.groups))
  temp <- sim.matrix
  temp <- temp[, order(temp[1,])]
  for (j in 2:(length(method2.groups) - 1)){
    sub.temp <- temp[, 1:(ncol(temp) - j + 1)]
    temp <- cbind(sub.temp[, order(sub.temp[j,])], temp[, (ncol(temp) - j + 2):ncol(temp), drop=FALSE])
  }
  sim.matrix <- temp
  melted.sim.matrix <- melt(sim.matrix)

  ggplot(data=melted.sim.matrix, aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill=value)) +
    theme_classic() +
    scale_fill_gradientn(colours=rev(grDevices::heat.colors(10))) +
    xlab(method1.name) +
    ylab(method2.name) +
    theme(plot.title=element_text(size=14, face='bold', hjust=0.5),
          legend.title=element_text(size=16),
          axis.title=element_text(size=16),
          axis.title.x=element_text(vjust=-0.5),
          axis.title.y=element_text(vjust=2.5),
          axis.text=element_text(size=10, face='bold'),
          axis.text.x=element_text(angle=45, hjust=1))
}
