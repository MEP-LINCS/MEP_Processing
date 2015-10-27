library("rmarkdown")

for(cellLine in c("MCF7")){
  for(ss in c("SS1")){
    render("MEP-LINCS_Analysis.Rmd", output_file = paste0("Mep-LINCS_Analysis_",cellLine,"_",ss,".html"), output_format = "html_document") 
  }
}
