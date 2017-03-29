#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"

library(synapseClient)
library(MEMA)
library(parallel)
library(stringr)
suppressPackageStartupMessages(library(dplyr))
library(optparse)

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-e", "--excelMetadata"), action="store_true", default=FALSE,
                help="Get metadata from Excel files instead of from the An! database"),
    make_option(c("-l", "--local"), action="store_true", default=FALSE,
                help="Use a local server instead of Synpase for file accesses"),
    make_option(c("-r", "--rawDataVersion"), type="character", default="v2",
                help="Raw data version from local server [default \"%default\"]")
  )
  parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
###Debug
cl <-list(options=list(verbose=TRUE,
                       excelMetadata=FALSE,
                       local=FALSE,
                       rawDataVersion="v2"),
          args=c("/lincs/share/lincs_user/LI8X00641",
                 "/lincs/share/lincs_user/LI8X00641/Analysis/LI8X00641_Level1.tsv"))
####
cl <- getCommandLineArgs()

barcodePath <- cl$args[1]
barcode <- gsub(".*/", "", barcodePath)
path <- gsub(barcode, "", barcodePath)
ofname <- cl$args[2]

opt <- cl$options
verbose <- opt$verbose
useAnnotMetadata <- !opt$excelMetadata
useSynapse <- !opt$local
rawDataVersion <- opt$rawDataVersion

if(useSynapse) synapseLogin()

scriptStartTime <- Sys.time()

if (verbose) message(paste("Processing plate:", barcode, "at", path, "\n"))

if (useAnnotMetadata) {
  #Build metadata file name list
  if(useSynapse){
    metadataq <- sprintf("select id from syn8466225 WHERE DataType='Metadata' AND Barcode='%s'",barcode)
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
  q <- sprintf("select id,name,Barcode,Level,Well,StainingSet,Location,Study from syn7800478 WHERE Level='0' AND Barcode='%s'",
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

if(verbose) message("Writing cell level data\n")
fwrite(cDT, file=ofname, sep = "\t", quote = FALSE)
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
                                     Level = "1")
  
  synFile <- synStore(synFile,
                      used=c(dataBWInfo$id, metadataTable$id),
                      forceVersion=FALSE)
}

if(verbose) message(paste("Elapsed time:", Sys.time()-scriptStartTime, "\n"))
