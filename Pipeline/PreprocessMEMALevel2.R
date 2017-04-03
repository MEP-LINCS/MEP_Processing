#!/usr/bin/env Rscript

#author: "Mark Dane"

library(MEMA)
library(parallel)
library(stringr)
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(optparse))

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

cl <-list(options=list(verbose=TRUE),
          args=c("LI8X00641",
                 "/tmp/LI8X00641_Level2.tsv"))
####
## cl <- getL2CommandLineArgs()

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
if (useSynapse) {
  synapseLogin()
  fileViewId <- 'syn7494072'
  level <- "1"
  levelQuery <- sprintf('SELECT id from %s WHERE Barcode="%s" AND Level="%s"',
                         fileViewId, barcode, level)
  levelRes <- synTableQuery(levelQuery)

  if (nrow(levelRes@values) > 1) {
    stop(sprintf("Found more than one Level 1 file for barcode %s", barcode))
  }
  
  dataPath <- getFileLocation(synGet(levelRes@values$id))

  imageIdQuery <- sprintf('SELECT id from %s WHERE Barcode="%s" AND DataType="ImageID"',
                        fileViewId, barcode)
  imageIdRes <- synTableQuery(imageIdQuery)
  
  if (nrow(imageIdRes@values) > 1) {
    stop(sprintf("Found more than one ImageID file for barcode %s", barcode))
  }
  
  imageIdPath <- getFileLocation(synGet(imageIdRes@values$id))
  
  } else {
  dataPath <- paste0(path,"/",barcode,"_Level1.tsv")
  imageIdPath <- paste0(path, "/",barcode, "_imageIDs.tsv")
}

cDT <- fread(dataPath)
omeroIds <- getOmeroIDs(imageIdPath)

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
spotDT <- merge(spotDT, omeroIds,
                by=c("WellIndex","ArrayRow","ArrayColumn"))

if(verbose) message("Writing spot level data\n")
writeTime<-Sys.time()
fwrite(data.table(spotDT), file=ofname, sep = "\t", quote=FALSE)

if(useSynapse){
  annotatedFolder <- synStore(Folder(name='Annotated', parentId="syn4215176"))
  synFile <- File(ofname, parentId=annotatedFolder@properties$id)
  synSetAnnotations(synFile) <- list(CellLine = unique(cDT$CellLine),
                                     Barcode = barcode,
                                     Study = unique(cDT$Study),
                                     Preprocess = "v1.8",
                                     DataType = "Quanititative Imaging",
                                     Consortia = "MEP-LINCS",
                                     Drug = unique(cDT$Drug),
                                     Segmentation = rawDataVersion,
                                     StainingSet = gsub("Layout.*","",unique(cDT$StainingSet)),
                                     Level = "2")
  
  synFile <- synStore(synFile,
                      used=c(levelRes@values$id, imageIdRes@values$id),
                      forceVersion=FALSE)
}

if(verbose) message(paste("Write time:", Sys.time()-writeTime,"\n"))
if(verbose) message(paste("Elapsed time:", Sys.time()-functionStartTime, "\n"))

