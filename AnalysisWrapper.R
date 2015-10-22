library("rmarkdown")

for(cellLine in c("PC3", "MCF7", "YAPC")){
  for(ss in c("SS1", "SS2", "SS3")){
    render("MEP-LINCS_Analysis.Rmd", output_file = paste0("Mep-LINCS_Analysis_",cellLine,"_",ss,".html"), output_format = "html_document") 
  }
}
