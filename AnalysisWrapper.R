library("rmarkdown")

PC3df <- data.frame(filePath="",
                    cellLine=rep(c("PC3"), 4),
                    ss=c("SS1", "SS2","SS3","SS2noH3"),
                    drug=c("none"),
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1", "v1"),
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(filePath="",
                     cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion=c("av1.6"),
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(filePath="",
                     cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(filePath="",
                       cellLine="MCF10A",
                       ss=c("SS1","SS2","SS3"),
                       drug=c("none"),
                       analysisVersion="av1.6",
                       rawDataVersion="v2",
                       stringsAsFactors=FALSE)

HCC1954Lapatinibdf <- data.frame(filePath="~/Documents/ME Watson/Lapatinib MEMAs/HCC1954/Lapatinib/SS6/",
                                 cellLine="HCC1954",
                                 ss=c("SS6"),
                                 drug=c("Lapatinib"),
                                 analysisVersion="av1.6",
                                 rawDataVersion="v2",
                                 stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf,HCC1954Lapatinibdf)

renderAnalysisReports <- function(x){
  filePath<- x[["filePath"]] 
  cellLine <- x[["cellLine"]]
  ss <- x[["ss"]]
  drug <- x[["drug"]]
  rawDataVersion <- x[["rawDataVersion"]]
  analysisVersion <- x[["analysisVersion"]]
  render("MEP-LINCS_Analysis.Rmd",
         output_file = paste0("./AnalysisReports/MEP-LINCS_Analysis_",
                              x[["cellLine"]],"_",
                              x[["ss"]],"_",
                              x[["rawDataVersion"]],"_",
                              x[["analysisVersion"]],".html"),
         output_format = "html_document")
}
tmp <- apply(ssDatasets[c(14),], 1, renderAnalysisReports)
