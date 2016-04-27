library("rmarkdown")

PC3df <- data.frame(cellLine=rep(c("PC3"), 3),
                    ss=c("SS1", "SS2","SS3"),
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1"),
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     analysisVersion=c("av1.6"),
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine=rep(c("MCF10A"), 3),
                       ss=c("SS1","SS2","SS3"),
                       analysisVersion="av1.6",
                       rawDataVersion=c("v2","v2","v2"),
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf)

x <- ssDatasets[c(1:3),]

renderCommonSignalNormEvalReports <- function(x){
  cellLine <- unique(x[["cellLine"]])
  k=128L
  verbose<-FALSE
  analysisVersion <- unique(x[["analysisVersion"]])
  render("MEP-LINCS_CommonSigNormEvalx2.Rmd",
         output_file = paste0("./NormEvalReports/MEP-LINCS_CommonSigNormEvalx2",
                              cellLine,"_",
                              analysisVersion,"_",
                              "k",k,".html"),
         output_format = "html_document")
}

tmp <- renderCommonSignalNormEvalReports(x=x)
