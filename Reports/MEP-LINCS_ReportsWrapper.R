library(knitr)
library(rmarkdown)

renderQACellReport <- function(studyName, path){
  render("MEP-LINCS_QACellLevel.Rmd",
         output_file = paste0(path,"/study/",studyName,"/Reports/MEP-LINCS_QA_Cell_",studyName,".html"),
         output_format = "html_document")
}

renderQASpotMEPReport <- function(studyName, path){
  render("MEP-LINCS_QASpotMEPLevel.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/MEP-LINCS_QA_SpotMEP_",studyName,".html"),
         output_format = "html_document")
}

renderAnalysisReport <- function(studyName, path){
  render("MEP-LINCS_Analysis.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/MEP-LINCS_Analysis_",studyName,".html"),
         output_format = "html_document")
}

renderSSCAnalysisReports <- function(studyName, path, fileViewSynID){
  render("MEP-LINCS_AnalysisSB_SSC.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/MEP-LINCS_AnalysisSB_",
                              studyName,".html"),
         output_format = NULL)
}

renderComparePipelines <- function(studyName,path){
  render("../archive/ComparePipelines.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/ComparePipelines_",studyName,".html"),
         output_format = NULL)
}

path <- "/lincs/share/lincs_user"
studyNames <- c("HMEC122L_SS1","HMEC122L_SS4","HMEC240L_SS1","HMEC240L_SS4","MCF10A_SS1","MCF10A_SS2","MCF10A_SS3")
studyNames <- c("HMEC122L_SS1","HMEC122L_SS4","HMEC240L_SS1","HMEC240L_SS4")
studyNames <- c("MCF10A_SS1","MCF10A_SS2","MCF10A_SS3")
studyNames <- c("au565_dmso_ss100","au565_lapatinib_ss100","hcc1954_dmso_ss100","hcc1954_lapatinib_ss100")

res <- lapply(studyNames, renderQACellReport, path=path)

path <- "/lincs/share/lincs_user/study"
res <- lapply(studyNames,renderQASpotMEPReport, path=path)

res <- lapply(studyNames,renderAnalysisReport, path=path)

# SSCStudyNames <- c("au565_dms0_ssc")
# #res <- renderSSCAnalysisReports("HMEC240L_SSC", path, fileViewSynID="syn7494072")
# res <- lapply(SSCStudyNames,renderSSCAnalysisReports, path=path, fileViewSynID="syn9612057")
# #res <- lapply(tolower(studyNames),renderComparePipelines, path=path)
