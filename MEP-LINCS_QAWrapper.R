library("rmarkdown")

analysisVersion <- "v1.4"
rawDataVersion <- "v1"
for(cellLine in c("MCF7","PC3","YAPC")[c(2)]){
  for(ss in c("SS1","SS2noH3","SS2","SS3","SS0")[c(2)]){
    render("MEP-LINCS_QA.Rmd", output_file = paste0("./QAReports/MEP-LINCS_QA_",cellLine,"_",ss,"_",rawDataVersion,"_",analysisVersion,".html"), output_format = "html_document") 
  }
}

