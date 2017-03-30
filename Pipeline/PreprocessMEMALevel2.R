#!/usr/bin/env Rscript

#author: "Mark Dane"

# Get the command line arguments and options
# returns a list with options and args elements
getL2CommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-l", "--local"), type="character", default=NULL,
                help="Path to local input data directory if not using Synpase.")
  )
  parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(stringr)
suppressPackageStartupMessages(library(optparse))

cl <-list(options=list(verbose=TRUE,
                       local="/lincs/share/lincs_user/LI8X00641/Analysis"),
          args=c("LI8X00641",
                 "/lincs/share/lincs_user/LI8X00641/Analysis/LI8X00641_Level2.tsv"))
####
cl <- getL2CommandLineArgs()

barcode <- cl$args[1]
ofname <- cl$args[2]

opt <- cl$options
verbose <- opt$verbose
if(is.null(opt$local)){
  useSynapse <- TRUE
} else {
  useSynapse <- FALSE
  path <- opt$local
}

if (verbose) message(paste("Summarizing cell to spot data for plate",barcode,"\n"))
functionStartTime<- Sys.time()
startTime<- Sys.time()
seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")

#Read in the plate's cell level data
if(useSynapse){
  stop("Synapse not supported in level 2 yet")
} else {
  cDT <- fread(paste0(path,"/",barcode,"_Level1.tsv"))
}

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
spotDT <- merge(spotDT,getOmeroIDs(path, barcode),by=c("WellIndex","ArrayRow","ArrayColumn"))

if(verbose) message("Writing spot level data\n")
writeTime<-Sys.time()
fwrite(data.table(spotDT), file=ofname, sep = "\t", quote=FALSE)

if(useSynapse){
  stop("Synapse not supported in level 2 yet")
}

if(verbose) message(paste("Write time:", Sys.time()-writeTime,"\n"))
if(verbose) message(paste("Elapsed time:", Sys.time()-functionStartTime, "\n"))

