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
    make_option(c("-i", "--inputPath"), type="character", default=NULL, metavar="PATH",
                help="Path to local input data directory or Synapse ID for a File View.")
  )
  parser <- OptionParser(usage = "%prog [options] STUDY OUTPUTFILE", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
###Debug
cl <-list(options=list(verbose=TRUE,
                       inputPath="/lincs/share/lincs_user/study/HMEC122L_SS1/Annotated/HMEC122L_SS1_Level3.tsv"),
          args=c("HMEC122L_SS1",
                 "/lincs/share/lincs_user/study/HMEC122L_SS1/Annotated/HMEC122L_SS1_Level4.tsv")
)
####
cl <- getL4CommandLineArgs()

studyName <- cl$args[1]
ofname <- cl$args[2]

opt <- cl$options
verbose <- opt$verbose
if(is.null(opt$local)){
  useSynapse <- TRUE
} else {
  useSynapse <- FALSE
  ifname <- opt$local
}

startTime <- Sys.time()
if(useSynapse){
  stop("Synapse not supported for level 4 data yet")
} else {
  l3DT <- fread(ifname)
}
#Summarize to the MEP_Drug level
mepDT <- preprocessLevel4(l3DT,seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin"))
#Add in the barcodes for each MEP_Drug
mepDT <- addBarcodes(dt3 = l3DT, dt4 = mepDT)
# Add a QA flag for spots with few replicates
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < 3

fwrite(mepDT, file=ofname, sep = "\t", quote=FALSE)
if(verbose) message("Writing level 4 file to disk\n")
if(useSynapse){
  stop("Synpase not supported in level 4 data yet")
}
if(verbose) message(paste("Elapsed time for ",studyName, "is", Sys.time()-startTime, "\n"))

