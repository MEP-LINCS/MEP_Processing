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

renderSSCAnalysisReports <- function(studyName, path){
  render("MEP-LINCS_AnalysisSB_SSC.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/MEP-LINCS_AnalysisSB_",
                              studyName,".html"),
         output_format = NULL)
}

studyNames <- c("HMEC122L_SS1","HMEC122L_SS4","HMEC240L_SS1","HMEC240L_SS4","MCF10A_SS1","MCF10A_SS2","MCF10A_SS3","au565_dmso_ss100","au565_lapatinib_ss100","hcc1954_dmso_ss100","hcc1954_lapatinib_ss100", "mcf10a_em_ssi")

path <- "/lincs/share/lincs_user"
res <- lapply(studyNames, renderQACellReport, path=path)

path <- "/lincs/share/lincs_user/study"
res <- lapply(studyNames,renderQASpotMEPReport, path=path)
res <- lapply(studyNames[8:12],renderAnalysisReport, path=path)

SSCStudyNames <- c("MCF10A_SSC","HMEC122L_SSC","HMEC240L_SSC")
#res <- lapply(SSCStudyNames,renderSSCAnalysisReports, path=path)
