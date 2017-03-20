library("rmarkdown")

renderAnalysisReports <- function(studyName, path){
  render("MEP-LINCS_AnalysisSB_SSC.Rmd",
         output_file = paste0(path,"/",studyName,"/Reports/MEP-LINCS_AnalysisSB_",
                              studyName,".html"),
         output_format = NULL)
}

r <- renderAnalysisReports(studyName="MCF10A_SSC", path="/lincs/share/lincs_user/study")
r <- renderAnalysisReports(studyName="HMEC122L_SSC", path="/lincs/share/lincs_user/study")
r <- renderAnalysisReports(studyName = "HMEC240L_SSC", path="/lincs/share/lincs_user/study")


