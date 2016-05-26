library("rmarkdown")

analysisVersion <- "v1.3"
for(cellLine in c("MCF7","PC3","YAPC")[2]){
  for(ss in c("SS1","SS2noH3","SS2","SS3","SS0")[5]){
    render("MEP-LINCS_QA_CellCount.Rmd", output_file = paste0("./QAReports/MEP-LINCS_QA_",cellLine,"_",analysisVersion,"_",ss,".html"), output_format = "html_document") 
  }
}

