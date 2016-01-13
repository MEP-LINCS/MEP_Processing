library("rmarkdown")
# 
analysisVersion <- "v1.3"
for(cellLine in c("MCF7","PC3","YAPC")){
  for(ss in c("SS1","SS2","SS3")){
    render("MEP-LINCS_Analysis.Rmd", output_file = paste0("./AnalysisReports/MEP-LINCS_Analysis_",cellLine,"_",analysisVersion,"_",ss,".html"), output_format = "html_document") 
  }
}

