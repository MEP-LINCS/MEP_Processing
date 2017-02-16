library("rmarkdown")

renderQACellReport <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_QACellLevel.Rmd",
         output_file = paste0(path,"/study/",datasetName,"/Reports/MEP-LINCS_QA_Cell_",datasetName,".html"),
         output_format = "html_document")
}

renderQASpotMEPReport <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_QASpotMEPLevel.Rmd",
         output_file = paste0(path,"/",datasetName,"/Reports/MEP-LINCS_QA_SpotMEP_",datasetName,".html"),
         output_format = "html_document")
}

renderAnalysisReport <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_Analysis.Rmd",
         output_file = paste0(path,"/",datasetName,"/Reports/MEP-LINCS_Analysis_",datasetName,".html"),
         output_format = "html_document")
}

path <- "/lincs/share/lincs_user"
datasetName="Longitudinal_study"
tmp <- renderQACellReport(path, datasetName)
path <- "/lincs/share/lincs_user/study"
tmp <- renderQASpotMEPReport(path, datasetName)
tmp <- renderAnalysisReport(path, datasetName)
