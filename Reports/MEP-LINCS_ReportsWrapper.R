library(knitr)
library(rmarkdown)
library(optparse)

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  parser <- OptionParser(usage = "%prog [options] PATH STUDYNAME")
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
if(!interactive()){
  cl <- getCommandLineArgs()
  path <- cl$args[1]
  studyName <- cl$args[2]
  
} else {
  path <- "/lincs/share/lincs_user"
  studyName <- "cama1_highserum_vehicle_mema"
}

render("MEP-LINCS_QACellLevel.Rmd",
       output_file = paste0(path,"/study/", studyName,"/Reports/MEP-LINCS_QA_Cell_",studyName,".html"),
       output_format = "html_document")

render("MEP-LINCS_QASpotMEPLevel.Rmd",
       output_file = paste0(path,"/study/",studyName,"/Reports/MEP-LINCS_QA_SpotMEP_",studyName,".html"),
       output_format = "html_document")

render("MEP-LINCS_Analysis.Rmd",
       output_file = paste0(path,"/study/",studyName,"/Reports/MEP-LINCS_Analysis_",studyName,".html"),
       output_format = "html_document")

# render("MEP-LINCS_AnalysisSB_SSC.Rmd",
#        output_file = paste0(path,studyName,"/Reports/MEP-LINCS_AnalysisSB_",
#                             studyName,".html"),
#        output_format = NULL)

