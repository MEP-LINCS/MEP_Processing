library("rmarkdown")

PC3df <- data.frame(cellLine=rep(c("PC3"), 3),
                    ss=c("SS1", "SS2","SS3"),
                    analysisVersion="av1.4",
                    rawDataVersion=c("v2","v2.1","v2.1"),
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     analysisVersion=c("av1.4", "av1.4","av1.4"),
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     analysisVersion="av1.4",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine=rep(c("MCF10A"), 2),
                       ss=c("SS1","SS3"),
                       analysisVersion="av1.4",
                       rawDataVersion=c("v2","v2"),
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf)

x <- ssDatasets[c(1:3),]
renderCommonSignalNormEvalReports <- function(x){
  cellLine <- unique(x[["cellLine"]])
  k=10
  verbose<-FALSE
  analysisVersion <- unique(x[["analysisVersion"]])
  render("MEP-LINCS_CommonSigNormEval.Rmd",
         output_file = paste0("./NormEvalReports/MEP-LINCS_CommonSigNormEval_",
                              cellLine,"_",
                              analysisVersion,".html"),
         output_format = "html_document")
}

tmp <- renderCommonSignalNormEvalReports(x=ssDatasets[1:3,])
