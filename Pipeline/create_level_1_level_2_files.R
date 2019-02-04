library(tidyverse)
library(optparse)

#Segment all images in a plate directory
getCommandLineArgs <- function(){
  parser <- OptionParser(usage = "%prog [options] INPUT_PATH OUTPUT_PATH PLATE_ID")
  arguments <- parse_args(parser, positional_arguments = 3)
}

#Specify the command line options
if(!interactive()){
  cl <- getCommandLineArgs()
  input_path <- cl$args[1]
  output_path <- cl$args[2]
  plate_ID <- cl$args[3]
} else {
  input_path <- "/lincs/share/lincs_user/"
  output_path <- "/lincs/share/lincs_user/"
  plate_ID <- "LI8V01171"
}

wells <- dir(paste0(input_path,plate_ID), pattern = "Well_._images")%>%
  str_remove(., "_images")
#wells <- wells[c(1,3)]
res <- lapply(wells, function(well){
  job_name <- paste0("l1_",str_extract(plate_ID,".{3}$"),str_extract(well,".{2}$"))
  sys_command <- paste0("srun --nodes 1 --cpus-per-task 3 --mem 30G --exclude=eppec-node6,eppec-node7 -o OutputCreateLevel1.txt -e ErrorCreateLevel1.txt --job-name=", job_name, "Rscript Preprocess_GT_Features_plate.R ",input_path, plate_ID)
  system(sys_command, wait = FALSE)
})
