library("rmarkdown")

analysisVersion <- "v1.3"
rawDataVersion <- "v1"
for(cellLine in c("MCF7","PC3","YAPC")){
  for(ss in c("SS1","SS2noH3","SS2","SS3","SS0")[c(1)]){
    render("MEP-LINCS_Analysis.Rmd", output_file = paste0("./AnalysisReports/MEP-LINCS_Analysis_",cellLine,"_",ss, "_",rawDataVersion,"_",analysisVersion,".html"), output_format = "html_document") 
  }
}

