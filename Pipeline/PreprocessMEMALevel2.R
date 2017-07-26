#!/usr/bin/env Rscript

#author: "Mark Dane"

library(MEMA)
library(parallel)
library(stringr)
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(optparse))
library(githubr)

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-i", "--inputPath"), type="character", default=NULL, metavar="PATH",
                help="Path to local input data directory or Synapse ID for a File View."),
    make_option(c("--synapseStore"), type="character", default=NULL, metavar="SYNAPSEID",
                help="Store output file in Synapse directory (provide Synapse ID of Folder to store).")
  )
  parser <- OptionParser(usage = "%prog [options] BARCODE OUTPUTFILE", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

cl <- getCommandLineArgs()
barcode <- cl$args[1]
ofname <- cl$args[2]
opt <- cl$options
verbose <- opt$verbose
if(file.exists(opt$inputPath)){
  useSynapse <- FALSE
} else {
  useSynapse <- TRUE
}

if (verbose) message(sprintf("Summarizing cell to spot data for plate %s", barcode))
functionStartTime<- Sys.time()
startTime<- Sys.time()
seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")

#Read in the plate's cell level data
if (useSynapse) {
  suppressMessages(synapseLogin())
  level <- "1"
  levelQuery <- sprintf('SELECT id,Segmentation,Preprocess,DataType,Study,Consortia,StainingSet,CellLine,Drug from %s WHERE Barcode="%s" AND Level="%s"',
                        opt$inputPath, barcode, level)
  levelRes <- synTableQuery(levelQuery)
  
  if (nrow(levelRes@values) > 1) {
    stop(sprintf("Found more than one Level 1 file for barcode %s", barcode))
  }
  
  dataPath <- getFileLocation(synGet(levelRes@values$id))
  
  imageIdQuery <- sprintf("SELECT id from %s WHERE Barcode='%s' AND DataType='ImageID'",
                          opt$inputPath, barcode)
  imageIdRes <- synTableQuery(imageIdQuery)
  
  clarionIdQuery <- sprintf('SELECT id from %s WHERE Barcode="%s" AND DataType="ClarionID"',
                            opt$inputPath, barcode)
  clarionIdRes <- synTableQuery(clarionIdQuery)
  
  if (nrow(imageIdRes@values) == 1) imageIdPath <- getFileLocation(synGet(imageIdRes@values$id))

  if (nrow(clarionIdRes@values) == 1)   clarionIdPath <- getFileLocation(synGet(clarionIdRes@values$id))
  
  
} else {
  dataPath <- paste0(opt$inputPath, "/",barcode, "_Level1.tsv")
  imageIdPath <- paste0(opt$inputPath, "/",barcode, "_imageIDs.tsv")
  clarionIdPath <- paste0(opt$inputPath, "/",barcode, "_clarionIDs.tsv")
}

cDT <- fread(dataPath)
if(exists("imageIdPath")){
if(file.exists(imageIdPath)) omeroIDs <- getOmeroIDs(imageIdPath)
}

if(exists("clarionIdPath")){
  if(file.exists(clarionIdPath)){
    clarionIDs <- getOmeroIDs(clarionIdPath)
    setnames(clarionIDs, "ImageID", "ClarionID")
    if (useSynapse) levelRes@values$Preprocess <- paste0(levelRes@values$Preprocess,".1")
  } 
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
if(sum(c("ArrayRow","ArrayColumn") %in% colnames(spotDT))==2) spotDT <- QASpotData(spotDT, lthresh = .6)

#Merge in Omero imageID links
if(exists("omeroIDs")) spotDT <- merge(spotDT, omeroIDs,
                                       by=c("WellIndex","ArrayRow","ArrayColumn"))
#Merge in Clarion imageID links
if(exists("clarionIDs")) spotDT <- merge(spotDT, clarionIDs,
                                         by=c("WellIndex","ArrayRow","ArrayColumn"))

if(verbose) message("Writing spot level data")
writeTime<-Sys.time()
#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}
for (j in colnames(spotDT)) data.table::set(spotDT, j = j, value = shorten(spotDT[[j]]))

fwrite(spotDT, file=ofname, sep = "\t", quote=FALSE)

if(!is.null(cl$options$synapseStore)){
  if(verbose) message(sprintf("Writing to Synapse Folder %s", opt$synapseStore))
  #get permlink from GitHub
  scriptLink <- "https://github.com/MEP-LINCS/MEP_Processing/"
  repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)
  if(!class(repo)=="try-error" ) scriptLink <- getPermlink(repo, "Pipeline/PreprocessMEMALevel2.R")
  synFile <- File(ofname, parentId=opt$synapseStore)
  synSetAnnotations(synFile) <- list(CellLine = levelRes@values$CellLine,
                                     Barcode = barcode,
                                     Study = levelRes@values$Study,
                                     Preprocess = levelRes@values$Preprocess,
                                     DataType = levelRes@values$DataType,
                                     Consortia = levelRes@values$Consortia,
                                     Drug = levelRes@values$Drug,
                                     Segmentation = levelRes@values$Segmentation,
                                     StainingSet = levelRes@values$StainingSet,
                                     Level = "2")
  
  synFile <- synStore(synFile,
                      used=c(levelRes@values$id, imageIdRes@values$id),
                      executed=scriptLink,
                      forceVersion=FALSE)
}

if(verbose) message(paste("Write time:", Sys.time()-writeTime,"\n"))
if(verbose) message(paste("Elapsed time:", Sys.time()-functionStartTime, "\n"))

