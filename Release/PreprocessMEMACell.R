#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
# 2/2017

library(synapseClient)
library(MEMA)
library(parallel)
library(stringr)
library(dplyr)

processCellCommandLine <- function(x, useAnnotMetadata=TRUE, rawDataVersion="v2", verbose="FALSE"){
  if(length(x)==0) stop("There must be a barcodePath argument in the command line call")
  barcodePath <- x[1]
  if((length(x)>1)) useAnnotMetadata <- x[2]
  if((length(x)>1)) rawDataVersion <- x[3]
  if((length(x)>1)) verbose <- x[4]
  list(barcodePath,useAnnotMetadata,rawDataVersion,verbose)
}

#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00528",TRUE, "v2", TRUE) #MCF10A_SS2
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00510",TRUE, "v2", TRUE) #MCF10A_SS3
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00751",TRUE,"v2",TRUE) #MCF10A_Neratinib_2
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00771",TRUE,"v2",TRUE)
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00765",TRUE,"v2",TRUE)
callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00641",TRUE,"v2",TRUE)
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00850",TRUE,"v2",TRUE)
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00850",FALSE,"v2",TRUE)
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI9V01610",FALSE,"v1",TRUE)
# callParams <- processCellCommandLine(commandArgs(trailingOnly = TRUE))
barcodePath <-callParams[[1]]
useAnnotMetadata <-as.logical(callParams[[2]])
rawDataVersion <-callParams[[3]]
verbose <- as.logical(callParams[[4]])

synapseLogin()

scriptStartTime<- Sys.time()

barcode <- gsub(".*/","",barcodePath)
path <- gsub(barcode,"",barcodePath)
if(verbose) message(paste("Processing plate:",barcode,"at",path,"\n"))
#Build metdata file name list
if(useAnnotMetadata){
  metadataFiles <- list(annotMetadata=paste0(path,barcode,"/Analysis/",barcode,"_an2omero.csv"))
} else {
  metadataFiles <- list(logMetadata = dir(paste0(path,barcode,"/Analysis"),pattern = "xml",full.names = TRUE),
                        spotMetadata = dir(paste0(barcodePath,"/Analysis"),pattern = "gal",full.names = TRUE),
                        wellMetadata =  dir(paste0(path,barcode,"/Analysis"),pattern = "xlsx",full.names = TRUE)
 )
}
#Get all metadata
metadata <- getMetadata(metadataFiles, useAnnotMetadata)

#Gather filenames of raw data
q <- sprintf("select id,name,Barcode,Level,Well,StainingSet,Location from syn7800478 WHERE Level='0' AND Barcode='%s'",
             barcode)
rawFiles <- synTableQuery(q)
dataBWInfo <- rawFiles@values

res <- lapply(dataBWInfo$id, synGet)

cellDataFilePaths <- unlist(lapply(res, getFileLocation))
# cellDataFilePaths <- dir(paste0(barcodePath,"/Analysis/",rawDataVersion), full.names = TRUE)
dataBWInfo$Path <- cellDataFilePaths

if(length(cellDataFilePaths)==0) stop("No raw data files found")

#Gather data from either CP or INCell
if("Nuclei" %in% dataBWInfo$Location) { #Cell Profiler data
  dtL <- getCPData(dataBWInfo = dataBWInfo, verbose=verbose)
} else if(any(grepl("96well",dataBWInfo$Path))) { #GE INCell data
  dtL <- getICData(cellDataFilePaths = cellDataFilePaths, endPoint488 = unique(metadata$Endpoint488),endPoint555 = unique(metadata$Endpoint555),endPoint647 = unique(metadata$Endpoint647), verbose=verbose)
} else {
  stop("Only data from Cell Profiler or INCell pipelines are supported")
}
#Clean up legacy issues in column names and some values
dtL <- lapply(dtL, cleanLegacyIssues)
#Filter our debris and cell clusters
dtL <- lapply(dtL, function(dt){
  filterObjects(dt,nuclearAreaThresh = 50, nuclearAreaHiThresh = 4000)})
#Add local XY and polar coordinates
dtL <- mclapply(dtL, addPolarCoords, mc.cores = detectCores())
#Add spot level normalizations for median intensities
dtL <- mclapply(dtL,spotNormIntensities, mc.cores = detectCores())
#Add adjacency values
dtL <- mclapply(dtL, calcAdjacency, mc.cores = detectCores())
#Merge the data with metadata
cDT <- merge(rbindlist(dtL),metadata,by=c("Barcode","Well","Spot"))
# Gate cells where possible
cDT <- gateCells(cDT)
#Write the annotated cell level files to disk
writeCellLevel(dt=cDT,path=path,barcode=barcode, verbose=verbose)
#Write the File Annotations for Synapse to tab-delimited file
write.table(c(
  CellLine = unique(cDT$CellLine),
  Preprocess = "v1.8",
  DataType = "Quanititative Imaging",
  Consortia = "MEP-LINCS",
  Drug = unique(cDT$Drug),
  Segmentation = rawDataVersion,
  StainingSet = gsub("Layout.*","",unique(cDT$StainingSet)),
  Level = "1"
),paste0(barcodePath, "/Analysis/", barcode,"_","Level1Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
if(verbose) message(paste("Elapsed time:", Sys.time()-scriptStartTime, "\n"))
