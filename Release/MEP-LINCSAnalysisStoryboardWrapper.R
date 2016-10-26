library("rmarkdown")



renderAnalysisReports <- function(datasetName){
  render("./MEP_LINCS/Release/MEP-LINCS_AnalysisSB_SSC.Rmd",
         output_file = paste0("../AnalysisReports/MEP-LINCS_AnalysisSB_",
                              datasetName,".html"),
         output_format = NULL)
}

r <- renderAnalysisReports(datasetName="MCF10A")
r <- renderAnalysisReports(datasetName="HMEC122L")
r <- renderAnalysisReports(datasetName="HMEC240L")


