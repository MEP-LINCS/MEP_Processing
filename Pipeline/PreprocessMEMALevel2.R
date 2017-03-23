#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2/1/2017

# Get the command line arguments and options
# returns a list with options and args elements
getL2CommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-l", "--local"), action="store_true", default=FALSE,
                help="Use a local server instead of Synpase for file accesses"),
    make_option(c("-w", "--writeFiles"), action="store_true", default=FALSE,
                help="Write output files to disk")
    )
  parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 1)
}

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(stringr)
suppressPackageStartupMessages(library(optparse))

cl <-list(options=list(verbose=TRUE,
                       local=FALSE,
                       writeFiles=TRUE),
          args="/lincs/share/lincs_user/LI8X00641")
####
cl <- getL2CommandLineArgs()

barcodePath <- cl$args
barcode <- gsub(".*/", "", barcodePath)
path <- gsub(barcode, "", barcodePath)

opt <- cl$options
verbose <- opt$verbose
useSynapse <- !opt$local
writeFiles <- opt$writeFiles

if (verbose) message(paste("Summarizing cell to spot data for plate",barcode,"at",barcodePath,"\n"))
functionStartTime<- Sys.time()
startTime<- Sys.time()
seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")

#Read in the plate's cell level data
if(useSynapse){
  stop("Synapse not supported in level 2 yet")
} else {
  cDT <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1.tsv"))
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

