#!/bin/bash Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
# 2/2017
#barcodePath <-commandArgs(trailingOnly = TRUE)
#barcodePath <- "/lincs/share/lincs_user/LI8X00771" #8 well An!, CP
#useAnnotMetadata=TRUE
#barcodePath <- "/lincs/share/lincs_user/LI8X00850" #8 well !An! CP
#useAnnotMetadata=FALSE
barcodePath <- "/lincs/share/lincs_user/lincs96well/LI9V01612" #96 well !An! IC
useAnnotMetadata <- FALSE
verbose <- TRUE
rawDataVersion <- "v1"

scriptStartTime<- Sys.time()
library(MEMA)#merge, annotate and normalize functions
library(data.table)#fast file reads, data merges and subsetting
library(parallel)#use multiple cores for faster processing
library(stringr)

barcode <- gsub(".*/","",barcodePath)
path <- gsub(barcode,"",barcodePath)
if (verbose) cat("Processing plate:",barcode,"at",path,"\n")
#Get all metadata
metadata <- getMetadata(barcode, path, useAnnotMetadata)
#Gather filenames of raw data
cellDataFilePaths <- dir(paste0(barcodePath,"/Analysis/",rawDataVersion), full.names = TRUE)
if(length(cellDataFilePaths)==0) stop("No raw data files found")
dataBWInfo <- data.table(Path=cellDataFilePaths,
                         Well=gsub("_","",str_extract(dir(paste0(barcodePath,"/Analysis/",rawDataVersion)),"_.*_")),
                         Location=str_extract(cellDataFilePaths,"Nuclei|Cytoplasm|Cells|Image"))
#Gather data from either CP or INCell
if("Nuclei" %in% dataBWInfo$Location) {
  dtL <- getCPData(dataBWInfo = dataBWInfo, verbose=verbose)
} else if(any(grepl("96well",dataBWInfo$Path))) {
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
dtL <- lapply(dtL, addPolarCoords)
#Add spot level normalizations for median intensities
dtL <- lapply(dtL,spotNormIntensities)
#Add adjacency values
dtL <- mclapply(dtL, calcAdjacency,mc.cores = detectCores())
#Merge the data with metadata
cDT <- merge(rbindlist(dtL),metadata,by=c("Barcode","Well","Spot"))
# Gate cells where possible
cDT <- gateCells(cDT)
#Write the annotated cell level files to disk
writeCellLevel(dt=cDT,path=path,barcode=barcode, verbose=TRUE)
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
cat("Elapsed time:", Sys.time()-scriptStartTime, "\n")



