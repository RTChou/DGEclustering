genes <- unique(read.table('geneids.txt', header=F, stringsAsFactors=F)$V1)

library(GOSemSim)
library(AnnotationDbi)
library(AnnotationHub)
hub <- AnnotationHub::.Hub("AnnotationHub",
        getAnnotationHubOption("URL"),
        './AnnotationHubCache',
        httr::use_proxy(Sys.getenv("http_proxy")),
        FALSE)
orgdb <- query(hub, 'Drosophila melanogaster')[['AH57972']]

semData <- godata(OrgDb=orgdb, ont='BP', computeIC=FALSE)
Egenes <- mapIds(orgdb, keys=genes, column='ENTREZID', keytype='FLYBASE', multiVals='first')
Egenes <- Egenes[!is.na(Egenes)]
test <- mgeneSim(Egenes, semData, measure='Wang')
