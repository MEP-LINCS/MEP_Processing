#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
# 2/1/2017

processSpotCommandLine <- function(x, rawDataVersion="v2", verbose="FALSE"){
  if(length(x)==0) stop("There must be a barcodePath argument in the command line call")
  barcodePath <- x[1]
  if((length(x)>1)) rawDataVersion <- x[2]
  if((length(x)>1)) verbose <- x[3]
  list(barcodePath,rawDataVersion,verbose)
}

#Debug: To be deleted
#Issues: 
#callParams <- processSpotCommandLine("/lincs/share/lincs_user/LI8X00771","v2",TRUE)
#callParams <- processSpotCommandLine("/lincs/share/lincs_user/LI8X00655","v2",TRUE)
#callParams <- processSpotCommandLine("/lincs/share/lincs_user/LI8X00850","v2",TRUE)
#callParams <- processSpotCommandLine("/lincs/share/lincs_user/LI8X00850","v2",TRUE)
#callParams <- processSpotCommandLine("/lincs/share/lincs_user/lincs96well/LI9V01610","v1",TRUE)
callParams <- processSpotCommandLine(commandArgs(trailingOnly = TRUE))
barcodePath <-callParams[[1]]
rawDataVersion <-callParams[[2]]
verbose <- as.logical(callParams[[3]])

barcode <- gsub(".*/","",barcodePath)
path <- gsub(barcode,"",barcodePath)
if (verbose) cat("Summarizing cell to spot data for plate",barcode,"at",barcodePath,"\n")
functionStartTime<- Sys.time()
startTime<- Sys.time()
seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")

#library(limma)#read GAL file and strsplit2
library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(stringr)

#Read in the plate's cell level data and annotations
cDT <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1.tsv"))
annotations <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1Annotations.tsv"),header = FALSE)

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

if(verbose) cat("Writing spot level data to disk\n")
writeTime<-Sys.time()
fwrite(data.table(spotDT), paste0(barcodePath, "/Analysis/", barcode,"_","SpotLevel.tsv"), sep = "\t", quote=FALSE)
if(verbose) cat("Write time:", Sys.time()-writeTime,"\n")

#Write the File Annotations for Synapse to tab-delimited file
write.table(c(
  CellLine = annotations$V2[annotations$V1=="CellLine"],
  Preprocess = annotations$V2[annotations$V1=="Preprocess"],
  DataType = annotations$V2[annotations$V1=="DataType"],
  Consortia = annotations$V2[annotations$V1=="Consortia"],
  Drug = annotations$V2[annotations$V1=="Drug"],
  Segmentation = annotations$V2[annotations$V1=="Segmentation"],
  StainingSet = annotations$V2[annotations$V1=="StainingSet"],
  Level = "Spot"),
  paste0(barcodePath, "/Analysis/", barcode,"_","SpotLevelAnnotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
if(verbose) cat("Elapsed time:", Sys.time()-functionStartTime, "\n")

