#!/usr/bin/env Rscript

#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 2017

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(RUVnormalize)
library(ruv)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(synapseClient))

getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-k", "--k"), type="integer", default=256,
                help="Number of factors to use in RUV normalization [default %default]",
                metavar="number"),
    make_option(c("-i", "--inputSynID"), type="character", default=NULL, metavar="synID",
                help="Synapse ID for a File View with input data."),
    make_option(c("--synapseStore"), type="character", default=NULL, metavar="SYNAPSEID",
                help="Store output file in Synapse directory (provide Synapse ID of Folder to store).")
  )
  parser <- OptionParser(usage = "%prog [options] STUDY STUDY [STUDY]", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = c(2, 3))
}

#Get the command line arguments and options
cl <- getCommandLineArgs()
#Convert the command line arguments to a list of study names
studyNameList <- as.list(cl$args)
#Parse the command line options
opt <- cl$options
verbose <- opt$verbose
k <- opt$k

startTime <- Sys.time()
suppressMessages(synapseLogin())
#Get annotations from the level 3 files
annotations <- rbindlist(lapply(studyNameList, function(studyName){
  levelQuery <- sprintf('SELECT id,Segmentation,Preprocess,Study,DataType,Consortia,StainingSet,CellLine,Drug from %s WHERE Study="%s" AND Level="%s"', cl$options$inputSynID, studyName, 3)
  levelRes <- synTableQuery(levelQuery)
  levelRes@values
}))
#Check if any annotations contain more than one value 
if(nrow(unique(annotations[,.(Segmentation,Preprocess,DataType,Consortia,Drug)]))>1) stop("Cannot combine studies with mixed annotations")
#Read the level 2 files for the studies
level <- "2"
dataPathsL <- lapply(studyNameList, function(studyName){
  levelQuery <- sprintf('SELECT id,Segmentation,Preprocess,Study,DataType,Consortia,StainingSet,CellLine,Drug from %s WHERE Study="%s" AND Level="%s"',
                        cl$options$inputSynID, studyName, level)
  levelRes <- synTableQuery(levelQuery)
  dataPaths <- sapply(levelRes@values$id,
                      function(x) getFileLocation(synGet(x)))
})

#RUV and loess normalize the common DAPI signals
l3 <- preprocessCommonSignalsLevel3(paths = dataPathsL, k, verbose)

#Get the level 3 raw and normalized data
level <- "3"
dataPathsL3 <- lapply(studyNameList, function(studyName){
  levelQuery <- sprintf('SELECT id,Segmentation,Preprocess,Study,DataType,Consortia,StainingSet,CellLine,Drug from %s WHERE Study="%s" AND Level="%s"',
                        cl$options$inputSynID, studyName, level)
  levelRes <- synTableQuery(levelQuery)
  dataPaths <- sapply(levelRes@values$id,
                      function(x) getFileLocation(synGet(x)))
})
l3L <- lapply(dataPathsL3, function(dataPathL3){
  dt <- fread(unique(paste0(dataPathL3)), showProgress = FALSE)
  return(dt)
})

#Make a list of datatables containing the staining set specific columns and BW and Spot
sssL <- lapply(l3L, function(dt, l3){
  specificColumns <- c("BW","Spot",setdiff(colnames(dt),colnames(l3)))
  dt[,specificColumns, with=FALSE]
},l3=l3)
sssDT <- rbindlist(sssL, fill=TRUE)
#merge the common signals and the staining set specific signals
l3C <-sssDT[l3,on=c("BW","Spot")]

#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}
for (j in colnames(l3C)) data.table::set(l3C, j = j, value = shorten(l3C[[j]]))

#Write the normalized level3 data to disk
message("Writing level 3 combined data to disk\n")
studyNameSSC <- gsub("_.*","_ssc",studyNameList[[1]])
if(verbose) message(paste("Writing level 3 file to disk\n"))
ofname <- paste0("/tmp/",studyNameSSC,"_Level3.tsv")
fwrite(l3C, file = ofname, sep = "\t", quote=FALSE)

if(!is.null(cl$options$synapseStore)) {
  #get permlink from GitHub
  scriptLink <- "https://github.com/MEP-LINCS/MEP_Processing/"
  repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)
  if(!class(repo)=="try-error" ) scriptLink <- getPermlink(repo, "Pipeline/PreprocessCombinedStudies.R")
  if(verbose) message(sprintf("Writing to Synapse Folder %s", cl$options$synapseStore))
  synFile <- File(ofname, parentId=opt$synapseStore)
  synSetAnnotations(synFile) <- list(CellLine = unique(annotations$CellLine),
                                     Segmentation = unique(annotations$Segmentation),
                                     Preprocess = unique(annotations$Preprocess),
                                     DataType = unique(annotations$DataType),
                                     Consortia = unique(annotations$Consortia),
                                     Drug = unique(annotations$Drug),
                                     Study = studyNameSSC,
                                     StainingSet = "SSC",
                                     Level = "3")
  synFile <- synStore(synFile,
                      used=names(unlist(dataPathsL)),
                      executed=scriptLink,
                      forceVersion=FALSE)
}

#Summarize to the MEP_Drug_Drug1Conc level
mepDT <- preprocessLevel4(l3C[,grep("Endpoint|StainingSet|395nm|488nm|555nm|640nm|750nm|1An$|2An$|3An$|LigandSet",colnames(l3C),value=TRUE,invert=TRUE), with=FALSE],seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin"))

#Add in the barcodes for each MEP_Drug
mepDT <- addBarcodes(dt3 = l3C, dt4 = mepDT)

# Add a QA flag for spots with few replicates
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < 3

#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}
for (j in colnames(mepDT)) data.table::set(mepDT, j = j, value = shorten(mepDT[[j]]))

message("Writing level 4 file to disk\n")
ofname <- paste0("/tmp/",studyNameSSC,"_Level4.tsv")
fwrite(mepDT, file = ofname, sep = "\t", quote=FALSE)

if(!is.null(cl$options$synapseStore)) {
  if(verbose) message(sprintf("Writing to Synapse Folder %s", cl$options$synapseStore))
  synFile <- File(ofname, parentId=opt$synapseStore)
  synSetAnnotations(synFile) <- list(CellLine = unique(annotations$CellLine),
                                     Segmentation = unique(annotations$Segmentation),
                                     Preprocess = unique(annotations$Preprocess),
                                     DataType = unique(annotations$DataType),
                                     Consortia = unique(annotations$Consortia),
                                     Drug = unique(annotations$Drug),
                                     Study = studyNameSSC,
                                     StainingSet = "SSC",
                                     Level = "4")
  synFile <- synStore(synFile,
                      used=names(unlist(dataPathsL)),
                      executed=scriptLink,
                      forceVersion=FALSE)
}

message("Running time:",Sys.time()-startTime)

