library("rmarkdown")

renderQACellReports <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_QACellLevel.Rmd",
         output_file = paste0(path,"/",datasetName,"/Reports/MEP-LINCS_QA_Cell_",datasetName,".html"),
         output_format = "html_document")
}

renderQASpotMEPReports <- function(path, datasetName){
  render("MEP_LINCS/Release/MEP-LINCS_QASpotMEPLevel.Rmd",
         output_file = paste0(path,"/",datasetName,"/Reports/MEP-LINCS_QA_SpotMEP_",datasetName,".html"),
         output_format = "html_document")
}

path <- "/lincs/share/lincs_user"
studyName="HMEC122L_SS4"
tmp <- renderQACellReports(path, studyName)
tmp <- renderQASpotMEPReports(path, studyName)


