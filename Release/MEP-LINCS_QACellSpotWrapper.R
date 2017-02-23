library("rmarkdown")

PC3df <- data.frame(datasetName=c("PC3_SS1","PC3_SS2","PC3_SS3","PC3_SS2noH3"),
                    cellLine=rep(c("PC3"), 4),
                    ss=c("SS1", "SS2","SS3","SS2noH3"),
                    drug=c("none"),
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1", "v1"),
                    limitBarcodes=8,
                    k=7,
                    calcAdjacency=TRUE,
                    writeFiles = TRUE,
                    mergeOmeroIDs = TRUE,
                    useJSONMetadata=TRUE,
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(datasetName=c("MCF7_SS1","MCF7_SS2","MCF7_SS3"),
                     cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     limitBarcodes=8,
                     k=7,
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useJSONMetadata=TRUE,
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(datasetName=c("YAPC_SS1","YAPC_SS2","YAPC_SS3"),
                     cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     limitBarcodes=8,
                     k=7,
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useJSONMetadata=TRUE,
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(datasetName=c("MCF10A_SS1","MCF10A_SS2","MCF10A_SS3"),
                       cellLine="MCF10A",
                       ss=c("SS1","SS2","SS3"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8,8),
                       k=c(7,7,7),
                       calcAdjacency=TRUE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useJSONMetadata=TRUE,
                       stringsAsFactors=FALSE)

watsonMEMAs <- data.frame(datasetName=c("HCC1954_DMSO","HCC1954_Lapatinib","AU565_DMSO","AU565_Lapatinib"),
                          cellLine=c("HCC1954","HCC1954","AU565","AU565"),
                          ss=c("SS6"),
                          drug=c("DMSO","Lapatinib"),
                          analysisVersion="av1.6",
                          rawDataVersion="v2",
                          limitBarcodes=c(8,2,2,2),
                          k=c(7,1,1,1),
                          calcAdjacency=TRUE,
                          writeFiles = TRUE,
                          mergeOmeroIDs = TRUE,
                          useJSONMetadata=FALSE,
                          stringsAsFactors=FALSE)

qualPlates <- data.frame(datasetName="MCF10A_Qual",
                         cellLine=c("MCF10A"),
                         ss=c("SS0"),
                         drug=c("none"),
                         analysisVersion="av1.6",
                         rawDataVersion="v2",
                         limitBarcodes=c(4),
                         k=c(0),
                         calcAdjacency=FALSE,
                         writeFiles = TRUE,
                         mergeOmeroIDs = TRUE,
                         useJSONMetadata=FALSE,
                         stringsAsFactors=FALSE)

ctrlPlates <- data.frame(datasetName="MCF10A_Ctrl",
                         cellLine=c("MCF10A"),
                         ss=c("SS0"),
                         drug=c("none"),
                         analysisVersion="av1.6",
                         rawDataVersion="v2",
                         limitBarcodes=c(1),
                         k=c(0),
                         calcAdjacency=FALSE,
                         writeFiles = TRUE,
                         mergeOmeroIDs = TRUE,
                         useJSONMetadata=FALSE,
                         stringsAsFactors=FALSE)

HMEC240L <- data.frame(datasetName=c("HMEC240L_SS1","HMEC240L_SS4"),
                       cellLine=c("HMEC240L"),
                       ss=c("SS1","SS4"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8),
                       k=c(7,7),
                       calcAdjacency=TRUE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useJSONMetadata=TRUE,
                       stringsAsFactors=FALSE)

HMEC122L <- data.frame(datasetName=c("HMEC122L_SS1","HMEC122L_SS4"),
                       cellLine=c("HMEC122L"),
                       ss=c("SS1","SS4"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8),
                       k=c(7,7),
                       calcAdjacency=TRUE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useJSONMetadata=TRUE,
                       stringsAsFactors=FALSE)
ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf,watsonMEMAs,qualPlates, ctrlPlates, HMEC240L, HMEC122L)

tcDataSet <- data.frame(datasetName=c("MCF10A_TC"),
                        cellLine=c("MCF10A"),
                        ss=c("SS4"),
                        drug=c("none"),
                        analysisVersion="av1.7",
                        rawDataVersion="v2",
                        limitBarcodes=c(3),
                        k=c(64),
                        calcAdjacency=TRUE,
                        writeFiles = TRUE,
                        mergeOmeroIDs = TRUE,
                        useAnnotMetadata=FALSE,
                        stringsAsFactors=FALSE)

Bornstein <- data.frame(datasetName=c("BornsteinOSC","BornsteinCal27"),
                        cellLine=c("OSC","Cal27"),
                        ss=c("SSA"),
                        drug=c("radiation"),
                        analysisVersion="av1.7",
                        rawDataVersion="v2",
                        limitBarcodes=c(2),
                        k=c(64),
                        calcAdjacency=TRUE,
                        writeFiles = TRUE,
                        mergeOmeroIDs = TRUE,
                        useAnnotMetadata=FALSE,
                        stringsAsFactors=FALSE)

Baylor <- data.frame(datasetName=c("Baylor1", "Baylor2","Baylor3", "Baylor4","Baylor5", "Baylor6"),
                     cellLine=c("LM2", "LM2","LM2","SUM159","SUM159","SUM159"),
                     ss=c("SSD"),
                     drug=c("unknown"),
                     analysisVersion="av1.7",
                     rawDataVersion="v2",
                     limitBarcodes=c(2,2,1,2,2,1),
                     k=c(64),
                     calcAdjacency=FALSE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = FALSE,
                     useAnnotMetadata=FALSE,
                     stringsAsFactors=FALSE)

Vertex <- data.frame(datasetName=c("Vertex1", "Vertex2"),
                     cellLine=c("LCSC-311"),
                     ss=c("SSE"),
                     drug=c("none"),
                     analysisVersion="av1.7",
                     rawDataVersion="v2",
                     limitBarcodes=c(8),
                     k=c(64),
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useAnnotMetadata=FALSE,
                     stringsAsFactors=FALSE)

MLDDataSet <- data.frame(datasetName=c("MCF10ANeratinib2","MCF10ADMSO2","MCF10AVorinostat","MCF10ATrametinib"),
                         cellLine=c("MCF10A"),
                         ss=c("SSF"),
                         drug=c("Neratinib","DMSO","Vorinostat","Ttametinib"),
                         analysisVersion="av1.7",
                         rawDataVersion="v2",
                         limitBarcodes=c(8),
                         k=c(64),
                         calcAdjacency=TRUE,
                         writeFiles = TRUE,
                         mergeOmeroIDs = TRUE,
                         useAnnotMetadata=TRUE,
                         stringsAsFactors=FALSE)


validations <- data.frame(datasetName=c("MCF10AHighRep1","MCF10AHighRep3"),
                          cellLine=c("MCF10A"),
                          ss=c("SS4"),
                          drug=c("none"),
                          analysisVersion="av1.7",
                          rawDataVersion="v2",
                          limitBarcodes=c(4),
                          k=c(0),
                          calcAdjacency=TRUE,
                          writeFiles = TRUE,
                          mergeOmeroIDs = TRUE,
                          useAnnotMetadata=FALSE,
                          stringsAsFactors=FALSE)


renderQACellReports <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_QACellLevel.Rmd",
         output_file = paste0(path,"/",datasetName,"/Reports/MEP-LINCS_QA_Cell_",datasetName,".html"),
         output_format = "html_document")
}

renderQASpotMEPReports <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_QASpotMEPLevel.Rmd",
         output_file = paste0(path,"/",datasetName,"/Reports/MEP-LINCS_QA_SpotMEP_",datasetName,".html"),
         output_format = "html_document")
}

path <- "/lincs/share/lincs_user"
datasetName="MCF10A_DMSO_2"
tmp <- renderQACellReports(path, datasetName)
tmp <- renderQASpotMEPReports(path, datasetName)


