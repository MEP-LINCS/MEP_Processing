PC3df <- data.frame(cellLine=rep(c("PC3"), 4),
                 ss=c("SS1", "SS2","SS3","SS2noH3"),
                 analysisVersion=1.4,
                 rawDataVersion=c("v2","v2.1","v2.1", "v1"),
                 limitBarcodes=8,
                 k=7,
                 calcAdjacency=TRUE,
                 writeFiles = TRUE,
                 mergeOmeroIDs = TRUE,
                 stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                    ss=c("SS1", "SS2","SS3"),
                    analysisVersion=1.4,
                    rawDataVersion=c("v2","v2","v2"),
                    limitBarcodes=8,
                    k=7,
                    calcAdjacency=TRUE,
                    writeFiles = TRUE,
                    mergeOmeroIDs = TRUE,
                    stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     analysisVersion=1.4,
                     rawDataVersion=c("v2","v2","v2"),
                     limitBarcodes=8,
                     k=7,
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine=rep(c("MCF10A"), 2),
                     ss=c("SS1","SS3"),
                     analysisVersion=1.4,
                     rawDataVersion=c("v2","v2"),
                     limitBarcodes=c(8,5),
                     k=c(7,4),
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     stringsAsFactors=FALSE)

datasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf)
k=7, limitBarcodes=8, analysisVersion="v1.4", rawDataVersion="v2", calcAdjacency=TRUE, writeFiles = TRUE, mergeOmeroIDs = TRUE, seNames=c("DN2N","SpotCellCount","EdU","MitoTracker","KRT","Fibrillarin")
