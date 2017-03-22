#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2/1/2017

#' Get the command line arguments and options
#' @return A list with options and args elements
#'@export
getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
                help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false",
                dest="verbose", help="Print little output"),
    make_option(c("-a", "--AnnotMetadata"), action="store_true", default=TRUE,
                help="use metadata from the An! database  [default]"),
    make_option(c("-e", "--ExcelMetadata"), action="store_false",
                dest="AnnotMetadata", help="use metadata from excel files"),
    make_option(c("-s", "--Synapse"), action="store_true", default=TRUE,
                help="use Synpase for file accesses  [default]"),
    make_option(c("-l", "--localServer"), action="store_false",
                dest="Synapse", help="use a local server for file access"),
    make_option(c("-w", "--writeFiles"), action="store_true", default=TRUE,
                help="write output files to disk  [default]"),
    make_option(c("-n", "--noWriteFiles"), action="store_false",
                dest="writeFiles", help="do not write output to disk"),
    make_option(c("-r", "--RawDataVersion"), type="character", default="v2",
                help="Raw data version from local server [default \"%default\"]")
  )
  parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 1)
}

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(stringr)
suppressPackageStartupMessages(library(optparse))


#Specify the command line options
###Debug
cl <-list(options=list(AnnotMetadata=TRUE,
                       verbose=TRUE,
                       Synapse=FALSE,
                       rawDataVersion="v2",
                       writeFiles=TRUE),
          args="/lincs/share/lincs_user/LI8X00641")
####
cl <- getCommandLineArgs()

barcodePath <- cl$args
barcode <- gsub(".*/", "", barcodePath)
path <- gsub(barcode, "", barcodePath)

opt <- cl$options
verbose <- opt$verbose
useSynapse <- opt$Synapse
writeFiles <- opt$writeFiles

if (verbose) message(paste("Summarizing cell to spot data for plate",barcode,"at",barcodePath,"\n"))
functionStartTime<- Sys.time()
startTime<- Sys.time()
seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")


#Read in the plate's cell level data and annotations
if(useSynapse){
  stop("Synapse not supported in level 2 yet")
} else {
  cDT <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1.tsv"))
}
#annotations <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1Annotations.tsv"),header = FALSE)

#Count the cells at each spot at the cell level as needed by createl3
cDT <- cDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]

#Add proportions for signals with multivariate gating and non-conforming gate values
addSpotProportions(cDT)

#Calculate proportions for binary gated signals
gatedSignals <- grep("Proportion", grep("Positive|High",colnames(cDT), value=TRUE), value=TRUE, invert=TRUE)
if(length(gatedSignals)>0){
  proportions <- cDT[,lapply(.SD, calcProportion),by="Barcode,Well,Spot", .SDcols=gatedSignals]
  setnames(proportions,
           grep("Gated",colnames(proportions),value=TRUE),
           paste0(grep("Gated",colnames(proportions),value=TRUE),"Proportion"))
}

#median summarize the rest of the signals to the spot level
signals <- summarizeToSpot(cDT, seNames)

if(exists("proportions")) {
  spotDT <- merge(signals,proportions)
} else {
  spotDT <- signals
}

#Set a threshold for the loess well level QA Scores
spotDT <- QASpotData(spotDT, lthresh = .6)

#Merge in Omero imageID links
spotDT <- merge(spotDT,getOmeroIDs(barcodePath),by=c("WellIndex","ArrayRow","ArrayColumn"))

if(writeFiles){
  if(verbose) message("Writing spot level data to disk\n")
  writeTime<-Sys.time()
  if(useSynapse){
    stop("Synapse not supported in level 2 yet")
  } else {
    fwrite(data.table(spotDT), paste0(barcodePath, "/Analysis/", barcode,"_","Level2.tsv"), sep = "\t", quote=FALSE)
  }
 
  if(verbose) message(paste("Write time:", Sys.time()-writeTime,"\n"))
}


if(verbose) message(paste("Elapsed time:", Sys.time()-functionStartTime, "\n"))

