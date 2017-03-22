#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2/2017

library(MEMA)
suppressPackageStartupMessages(library(optparse))


#' Get the command line arguments and options
#' @return A list with options and args elements
#'@export
getL4CommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
                help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false",
                dest="verbose", help="Print little output"),
    make_option(c("-s", "--Synapse"), action="store_true", default=TRUE,
                help="use Synpase for file accesses  [default]"),
    make_option(c("-l", "--localServer"), action="store_false",
                dest="Synapse", help="use a local server for file access"),
    make_option(c("-w", "--writeFiles"), action="store_true", default=TRUE,
                help="write output files to disk  [default]"),
    make_option(c("-n", "--noWriteFiles"), action="store_false",
                dest="writeFiles", help="do not write output to disk")
  )
  parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
###Debug
cl <-list(options=list(verbose=TRUE,
                       Synapse=FALSE,
                       writeFiles=FALSE),
          args=c("HMEC122L_SS1", "/lincs/share/lincs_user/study"))
####

startTime <- Sys.time()
#Get and decode command line arguments and options
cl <- getL4CommandLineArgs()
studyName <- cl$args[1]
path <- cl$args[2]

opt <- cl$options
verbose <- opt$verbose
useSynapse <- opt$Synapse
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

