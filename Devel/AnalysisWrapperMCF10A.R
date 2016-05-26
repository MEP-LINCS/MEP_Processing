library("rmarkdown")

MCF10Adf <- data.frame(cellLine=rep(c("MCF10A"), 2),
                       ss=c("SS1","SS3"),
                       analysisVersion="av1.4",
                       rawDataVersion=c("v2","v2"),
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(MCF10Adf)

renderAnalysisReports <- function(x){
  cellLine <- x[["cellLine"]]
  ss <- x[["ss"]]
  rawDataVersion <- x[["rawDataVersion"]]
  analysisVersion <- x[["analysisVersion"]]
  render("MEP-LINCS_AnalysisMCF10A.Rmd",
         output_file = paste0("./AnalysisReports/MEP-LINCS_Analysis_",
                              x[["cellLine"]],"_",
                              x[["ss"]],"_",
                              x[["rawDataVersion"]],"_",
                              x[["analysisVersion"]],".html"),
         output_format = "html_document")
}

tmp <- apply(ssDatasets[c(2),], 1, renderAnalysisReports)
