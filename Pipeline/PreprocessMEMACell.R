#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
# 2/2017

library(synapseClient)
library(MEMA)
library(parallel)
library(stringr)
library(dplyr)

processCellCommandLine <- function(x, useAnnotMetadata=TRUE, useSynapse=TRUE,
                                   rawDataVersion="v2", verbose="FALSE"){
  if(length(x)==0) stop("There must be a barcodePath argument in the command line call")
  barcodePath <- x[1]
  if (length(x) > 1) useAnnotMetadata <- x[2]
  if (length(x) > 1) useSynapse <- x[3]
  if (length(x) > 1) rawDataVersion <- x[4]
  if (length(x) > 1) verbose <- x[5]
  
  list(barcodePath, useAnnotMetadata, useSynapse, rawDataVersion, verbose)
}

#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00528",TRUE, "v2", TRUE) #MCF10A_SS2
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00510",TRUE, "v2", TRUE) #MCF10A_SS3
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00751",TRUE,"v2",TRUE) #MCF10A_Neratinib_2
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00771",TRUE,"v2",TRUE)
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00765",TRUE,"v2",TRUE)
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00641", TRUE, TRUE, "v2", TRUE)
callParams <- processCellCommandLine("/lincs/share/lincs_user/LI8X00831",TRUE,FALSE,"v2",TRUE)
#callParams <- processCellCommandLine("/lincs/share/lincs_user/LI9V01610",FALSE,"v1",TRUE)
# callParams <- processCellCommandLine(commandArgs(trailingOnly = TRUE))

barcodePath <- callParams[[1]]
useAnnotMetadata <- as.logical(callParams[[2]])
useSynapse <- as.logical(callParams[[3]])
rawDataVersion <- callParams[[4]]
verbose <- as.logical(callParams[[5]])

if(useSynapse) synapseLogin()

scriptStartTime <- Sys.time()

barcode <- gsub(".*/", "", barcodePath)
path <- gsub(barcode, "", barcodePath)

if (verbose) message(paste("Processing plate:", barcode, "at", path, "\n"))

if (useAnnotMetadata) {
  #Build metadata file name list
  if(useSynapse){
    metadataq <- sprintf("select id from syn8466225 WHERE DataType='Metadata' AND Barcode='%s'",
                         barcode)
    metadataFiles <- list(annotMetadata=synTableQuery(metadataq)@values$id)
    metadataTable <- synTableQuery(metadataq)@values
    metadataFiles <- lapply(metadataTable$id, synGet)
    metadataFiles <- list(annotMetadata=getFileLocation(metadataFiles[[1]]))
  } else {
    metadataFiles <- list(annotMetadata=paste0(barcodePath,"/Analysis/",barcode,"_an2omero.csv"))
  }
  
} else {
  metadataFiles <- list(logMetadata = dir(paste0(path,barcode,"/Analysis"),pattern = "xml",full.names = TRUE),
                        spotMetadata = dir(paste0(barcodePath,"/Analysis"),pattern = "gal",full.names = TRUE),
                        wellMetadata =  dir(paste0(path,barcode,"/Analysis"),pattern = "xlsx",full.names = TRUE)
  )
}

#Get all metadata
metadata <- getMetadata(metadataFiles, useAnnotMetadata)

#Gather filenames and metadata of level 0 files
if(useSynapse){
  q <- sprintf("select id,name,Barcode,Level,Well,StainingSet,Location from syn7800478 WHERE Level='0' AND Barcode='%s'",
               barcode)
  rawFiles <- synTableQuery(q)
  dataBWInfo <- rawFiles@values
  
  # Download raw files, or get from cache if already downloaded
  res <- lapply(dataBWInfo$id, synGet)
  
  # Get the path on disk to each file
  cellDataFilePaths <- unlist(lapply(res, getFileLocation))
  dataBWInfo$Path <- cellDataFilePaths
} else {
  cellDataFilePaths <- dir(paste0(barcodePath,"/Analysis/",rawDataVersion), full.names = TRUE)
  if(length(cellDataFilePaths)==0) stop("No raw data files found")
  dataBWInfo <- data.table(Path=cellDataFilePaths,
                           Well=gsub("_","",str_extract(dir(paste0(barcodePath,"/Analysis/",rawDataVersion)),"_.*_")),
                           Location=str_extract(cellDataFilePaths,"Nuclei|Cytoplasm|Cells|Image"))
}

if(length(cellDataFilePaths) == 0) stop("No raw data files found")

#Gather data from either Cell Profiler or INCell
if("Nuclei" %in% dataBWInfo$Location) { #Cell Profiler data
  dtL <- getCPData(dataBWInfo = dataBWInfo, verbose=verbose)
} else if (any(grepl("96well", dataBWInfo$Path))) { #GE INCell data
  dtL <- getICData(cellDataFilePaths = cellDataFilePaths, 
                   endPoint488 = unique(metadata$Endpoint488),
                   endPoint555 = unique(metadata$Endpoint555),
                   endPoint647 = unique(metadata$Endpoint647),
                   verbose=verbose)
} else {
  stop("Only data from Cell Profiler or INCell pipelines are supported")
}

#Clean up legacy issues in column names and some values
dtL <- lapply(dtL, cleanLegacyIssues)

# Filter our debris and cell clusters
dtL <- lapply(dtL, function(dt){
  filterObjects(dt,nuclearAreaThresh = 50, nuclearAreaHiThresh = 4000)})

# Add local XY and polar coordinates
# dtL <- mclapply(dtL, addPolarCoords, mc.cores = detectCores())
dtL <- lapply(dtL, addPolarCoords)

# Add spot level normalizations for median intensities
# dtL <- mclapply(dtL,spotNormIntensities, mc.cores = detectCores())
dtL <- lapply(dtL,spotNormIntensities)

# Add adjacency values
# dtL <- mclapply(dtL, calcAdjacency, mc.cores = detectCores())
dtL <- lapply(dtL, calcAdjacency)

# Merge the data with metadata
cDT <- merge(rbindlist(dtL),metadata,by=c("Barcode","Well","Spot"))

# Gate cells where possible
cDT <- gateCells(cDT)

# Write the annotated cell level files to disk
ofname <- paste0(path, barcode,"/Analysis/",barcode, "_", "Level1.tsv")
if(useSynapse){
  annotatedFolder <- synStore(Folder(name='Annotated', parentId="syn8466337"))
  synFile <- File(ofname, parentId=annotatedFolder@properties$id)
  synSetAnnotations(synFile) <- list(CellLine = unique(cDT$CellLine),
                                     Preprocess = "v1.8",
                                     DataType = "Quanititative Imaging",
                                     Consortia = "MEP-LINCS",
                                     Drug = unique(cDT$Drug),
                                     Segmentation = rawDataVersion,
                                     StainingSet = gsub("Layout.*","",unique(cDT$StainingSet)),
                                     Level = "1")
  
  synFile <- synStore(synFile,
                      used=c(dataBWInfo$id, metadataTable$id),
                      forceVersion=FALSE)
  
} else {
  fwrite(cDT, file=ofname, sep = "\t", quote = FALSE)
}


# #Write the File Annotations for Synapse to tab-delimited file
# write.table(c(
#   CellLine = unique(cDT$CellLine),
#   Preprocess = "v1.8",
#   DataType = "Quanititative Imaging",
#   Consortia = "MEP-LINCS",
#   Drug = unique(cDT$Drug),
#   Segmentation = rawDataVersion,
#   StainingSet = gsub("Layout.*","",unique(cDT$StainingSet)),
#   Level = "1"
# ),paste0(barcodePath, "/Analysis/", barcode,"_","Level1Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)

if(verbose) message(paste("Elapsed time:", Sys.time()-scriptStartTime, "\n"))
