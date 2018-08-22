#' @export
DGE.clust <- function(expressions, annotations, cluster.method='intego', nb.group, genclust.priori=NULL, nb.generation=500, LIM.ASSO = 4, LIM.COR = 0.5){
  nb.dim.ex = nrow(expressions)
  nb.dim.an = min((nrow(annotations) - 1), (ncol(annotations) - 1))

  annotations.sep = Integration(annotations, expressions, nb.dim.ex, LIM.ASSO, LIM.COR)
  annotations.sep <- apply(annotations.sep, 2, as.factor)
  MCA = MCAsimple(annotations.sep)[, 1:nb.dim.an]

  if (cluster.method == 'intego'){
    DIST <- dist(MCA, diag=TRUE, upper=TRUE)
    groups <- clustering(DIST, mode='Classification', nb.group=nb.group)
    res <- list(groups, annotations.sep, MCA)
  }

  else {
    # input for gensclust does not need to be scaled (i.e. MCA is already scaled)
    # generate MCA input for Genclust (i.e. input is GO terms integrated)
    line1 <-
      data.frame(as.list(c(
        nrow(MCA), ncol(MCA), 10, rep(NA, ncol(MCA) - 2)
      )))
    colnames(line1) <- seq(1, ncol(line1))

    line2 <- list('gene')
    for (i in seq(1, ncol(MCA))) {
      line2[i + 1] <- paste('f', i, sep = '')
    }
    line2 <- data.frame(line2)
    colnames(line2) <- colnames(line1)

    gen.input <- NULL
    gen.input <- cbind(as.data.frame(rownames(MCA)), as.data.frame(MCA))
    colnames(gen.input) <- colnames(line1)

    package.dir <- system.file(package="DGEclustering")
    gen.input <- rbind(line1, line2, gen.input)
    write.table(
      gen.input,
      file = paste(packae.dir, 'src/genclust_sig_data.tsv', sep='/'),
      sep = '\t',
      col.names = FALSE,
      row.names = FALSE,
      na = '',
      quote = FALSE
    )

    # generate initialization file
    filepath <- paste(packae.dir, 'src/out.tmp', sep='/')
    fileConn <- file(filepath)
    write('#Comment\n#Comment',
          filepath,
          append = FALSE,
          sep = '\n')
    if (is.null(genclust.priori) || genclust.priori=='random'){
        set.seed(123)
        ran <- split(seq(1, nrow(MCA)), sample(1:nb.group, nrow(MCA), replace = TRUE))
        for (i in 1:length(ran)) {
          write(paste(length(ran[[i]]), paste(ran[[i]], collapse = ' '), sep = '\t'),
                filepath,
                append = TRUE,
                sep = '\n')
        }
    }
    else {
      DIST <- dist(MCA, diag = TRUE, upper = TRUE)
      groups <- clustering(DIST, mode='Classification', nb.group=nb.group)
      rownames(gen.input) <- seq(1, nrow(gen.input))
      for (i in 1:nb.group){
        temp <- c()
        for (j in 1:length(groups[[i]])){
          temp <- c(temp, as.numeric(rownames(gen.input)[gen.input[, 1] == groups[[i]][j]]) - 2)
        }
        write(paste(length(groups[[i]]), paste(temp, collapse = ' '), sep = '\t'),
              filepath,
              append = TRUE,
              sep = '\n')
      }
    }
    close(fileConn)

    # run genclust
    system(
      paste(paste(packae.dir, 'src/genclust_sig_data.tsv', sep='/'),
        'genclust_sig_data.tsv',
        nb.group,
        nb.generation,
        'genclust_out.txt',
        0,
        0
      ), ignore.stdout = TRUE
    )

    # import genclust result
    gen.out <-
      read.table(paste(packae.dir, 'src/genclust_out.txt', sep='/'),
                 header = FALSE,
                 sep = '\t')
    gen.out[] <- lapply(gen.out, as.character)

    x <- list()
    g <- 0
    for (row in 1:(nrow(gen.out) - 1)) {
      string.list <- unlist(strsplit(gen.out[row, 1], ' '))
      if (string.list[1] == 'CLUSTER') {
        if (g > 0) {
          x <- list.append(x, temp)
        }
        g <- g + 1
        temp <- c()
        i <- 1
      }
      else{
        temp[i] <- string.list[2]
        i <- i + 1
      }
    }
    groups <- list.append(x, temp)
    names(groups) = paste('Group', 1:g, sep = '.')
    res <- list(groups)
  }

  names(res) <- c('groups', 'annotations.sep', 'MCA')
  return(res)
}
