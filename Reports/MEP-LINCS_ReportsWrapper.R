library(knitr)
library(rmarkdown)

renderQACellReport <- function(path, studyName){
  render("MEP-LINCS_QACellLevel.Rmd",
         output_file = paste0(path,"/study/",studyName,"/Reports/MEP-LINCS_QA_Cell_",studyName,".html"),
         output_format = "html_document")
}

renderQASpotMEPReport <- function(path, studyName){
  render("MEP-LINCS_QASpotMEPLevel.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/MEP-LINCS_QA_SpotMEP_",studyName,".html"),
         output_format = "html_document")
}

renderAnalysisReport <- function(path, studyName){
  render("MEP-LINCS_Analysis.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/MEP-LINCS_Analysis_",studyName,".html"),
         output_format = "html_document")
}

path <- "/lincs/share/lincs_user"
studyName="HMEC122L_SS1"
#tmp <- renderQACellReport(path, studyName)
path <- "/lincs/share/lincs_user/study"
#tmp <- renderQASpotMEPReport(path, studyName)
tmp <- renderAnalysisReport(path, studyName)
