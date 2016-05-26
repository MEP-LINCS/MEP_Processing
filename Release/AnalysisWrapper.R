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
                       analysisVersion="av1.6",
                       rawDataVersion="v2",
                       stringsAsFactors=FALSE)

watsonMEMAs <- data.frame(cellLine=c("HCC1954","HCC1954","AU565","AU565"),
                          ss=c("SS6"),
                          drug=c("DMSO","Lapatinib"),
                          analysisVersion="av1.6",
                          rawDataVersion="v2",
                          stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf,watsonMEMAs)

renderAnalysisReports <- function(x){
  cellLine <- x[["cellLine"]]
  ss <- x[["ss"]]
  drug <- x[["drug"]]
  rawDataVersion <- x[["rawDataVersion"]]
  analysisVersion <- x[["analysisVersion"]]
  render("./MEP_LINCS/MEP-LINCS_Analysis.Rmd",
         output_file = paste0("../AnalysisReports/MEP-LINCS_Analysis_",
                              x[["cellLine"]],"_",
                              x[["drug"]],"_",
                              x[["ss"]],"_",
                              x[["rawDataVersion"]],"_",
                              x[["analysisVersion"]],".html"),
         output_format = "html_document")
}
tmp <- apply(ssDatasets[c(16),], 1, renderAnalysisReports)
