library("rmarkdown")



renderAnalysisReports <- function(cellLine){
  render("./MEP_LINCS/Release/MEP-LINCS_AnalysisSB_SSC.Rmd",
         output_file = paste0("../AnalysisReports/MEP-LINCS_AnalysisSB_",
                              cellLine,"_SSC_",".html"),
         output_format = NULL)
}

r <- renderAnalysisReports(cellLine="MCF10A")
#r <- renderAnalysisReports(cellLine="HMEC122L")
#r <- renderAnalysisReports(cellLine="HMEC240L")



