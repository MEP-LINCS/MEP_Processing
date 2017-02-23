library("rmarkdown")

PC3df <- data.frame(cellLine=rep(c("PC3"), 4),
                    ss=c("SS1", "SS2","SS3","SS2noH3"),
                    drug=c("none"),
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1", "v1"),
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine="MCF10A",
                       ss=c("SS1","SS2","SS3"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       stringsAsFactors=FALSE)

watsonMEMAs <- data.frame(cellLine=c("HCC1954","HCC1954","AU565","AU565"),
                          ss=c("SS6"),
                          drug=c("DMSO","Lapatinib"),
                          analysisVersion="av1.6",
                          rawDataVersion="v2",
                          stringsAsFactors=FALSE)


HMEC240L <- data.frame(cellLine=c("HMEC240L"),
                       ss=c("SS1","SS4"),
                       drug="none",
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       stringsAsFactors=FALSE)

HMEC122L <- data.frame(cellLine=c("HMEC122L"),
                       ss=c("SS1","SS4"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf,watsonMEMAs, HMEC240L, HMEC122L)

Bornstein <- data.frame(datasetName=c("BornsteinOSC","BornsteinCal27"),
                        cellLine=c("OSC","Cal27"),
                        ss=c("SSA"),
                        drug=c("radiation"),
                        analysisVersion="av1.7",
                        rawDataVersion="v2",
                        useAnnotMetadata=FALSE,
                        stringsAsFactors=FALSE)

Vertex <- data.frame(datasetName=c("Vertex1","Vertex2"),
                     cellLine=c("LCSC-311"),
                     ss=c("SSE"),
                     drug=c("none"),
                     analysisVersion="av1.7",
                     rawDataVersion="v2",
                     stringsAsFactors=FALSE)

Baylor <- data.frame(datasetName=c("Baylor1", "Baylor2","Baylor3", "Baylor4","Baylor5", "Baylor6"),
                     cellLine=c("LM2", "LM2","LM2","SUM159","SUM159","SUM159"),
                     ss=c("SSD"),
                     drug=c("unknown"),
                     analysisVersion="av1.7",
                     rawDataVersion="v2",
                     stringsAsFactors=FALSE)

MLDDataSet <- data.frame(datasetName=c("MCF10ANeratinib2","MCF10ADMSO2","MCF10AVorinostat","MCF10ATrametinib"),
                         cellLine=c("MCF10A"),
                         ss=c("SSF"),
                         drug=c("Neratinib","DMSO","Vorinostat","Trametinib"),
                         analysisVersion="av1.7",
                         rawDataVersion="v2",
                         stringsAsFactors=FALSE)

validations <- data.frame(datasetName=c("MCF10AHighRep1", "MCF10AHighRep3"),
                          cellLine=c("MCF10A"),
                          ss=c("SS4"),
                          drug=c("none"),
                          analysisVersion="av1.7",
                          rawDataVersion="v2",
                          stringsAsFactors=FALSE)


renderAnalysisReports <- function(x){
  datasetName <- x[["datasetName"]]
  cellLine <- x[["cellLine"]]
  ss <- x[["ss"]]
  drug <- x[["drug"]]
  rawDataVersion <- x[["rawDataVersion"]]
  analysisVersion <- x[["analysisVersion"]]
  render("./MEP_LINCS/Release/MEP-LINCS_Analysis.Rmd",
         output_file = paste0("../AnalysisReports/MEP-LINCS_Analysis_",
                              x[["datasetName"]],"_",x[["ss"]],"_",".html"),
         output_format = "html_document")
}

tmp <- apply(MLDDataSet[1:2,], 1, renderAnalysisReports)


