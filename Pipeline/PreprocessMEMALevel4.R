#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2/2017

library(MEMA)
suppressPackageStartupMessages(library(optparse))

# Get the command line arguments and options
# returns a list with options and args elements
getL4CommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-l", "--local"), action="store_true", default=FALSE,
                help="Use a local server instead of Synpase for file accesses"),
    make_option(c("-w", "--writeFiles"), action="store_true", default=FALSE,
                help="Write output files to disk")
  )
  parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
###Debug
cl <-list(options=list(verbose=TRUE,
                       local=TRUE,
                       writeFiles=FALSE),
          args=c("HMEC122L_SS1",  "/lincs/share/lincs_user/study")
)
####
cl <- getL4CommandLineArgs()

studyName <- cl$args[1]
path <- cl$args[2]

opt <- cl$options
verbose <- opt$verbose
useSynapse <- !opt$local
writeFiles <- opt$writeFiles

startTime <- Sys.time()
if(useSynapse){
  stop("Synapse not supported for level 4 data yet")
} else {
  l3DT <- fread(paste0(path,"/",studyName,"/Annotated/",studyName,"_Level3.tsv"))
}
#Summarize to the MEP_Drug level
mepDT <- preprocessLevel4(l3DT,seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin"))
#Add in the barcodes for each MEP_Drug
mepDT <- addBarcodes(dt3 = l3DT, dt4 = mepDT)
# Add a QA flag for spots with few replicates
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < 3

if(writeFiles){
  if(verbose) message("Writing level 4 file to disk\n")
  if(useSynapse){
    stop("Synpase not supported in level 4 data yet")
  } else
    fwrite(mepDT, paste0(path, "/",studyName, "/Annotated/", studyName,"_Level4.tsv"), sep = "\t", quote=FALSE)
}

if(verbose) message(paste("Elapsed time for ",studyName, "is", Sys.time()-startTime, "\n"))

