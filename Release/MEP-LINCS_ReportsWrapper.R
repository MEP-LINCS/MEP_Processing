library("rmarkdown")

renderQACellReport <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_QACellLevel.Rmd",
         output_file = paste0(path,"/",datasetName,"/Reports/MEP-LINCS_QA_Cell_",datasetName,".html"),
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
datasetName="MCF10A_Neratinib_2"
tmp <- renderQACellReport(path, datasetName)
tmp <- renderQASpotMEPReport(path, datasetName)
tmp <- renderAnalysisReport(path, datasetName)
