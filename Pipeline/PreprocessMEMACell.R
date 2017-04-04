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
    make_option(c("-i", "--inputPath"), type="character", default=NULL, metavar="PATH",
                help="Path to local input data directory or Synapse ID for a File View."),
    make_option(c("-r", "--rawDataVersion"), type="character", default="v2",
                help="Raw data version from local server [default \"%default\"]"),
    make_option(c("--synapseStore"), type="character", default=NULL, metavar="SYNAPSEID",
                help="Store output file in Synapse directory (provide Synapse ID of Folder to store).")
  )
  parser <- OptionParser(usage = "%prog [options] BARCODE OUTPUTFILE", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
###Debug
cl <-list(options=list(verbose=TRUE,
                       excelMetadata=FALSE,
                       local="/lincs/share/lincs_user/LI8X00641/Analysis",
                       rawDataVersion="v2"),
          args=c("LI8X00641",
                 "/lincs/share/lincs_user/LI8X00641/Analysis/LI8X00641_Level1.tsv"))
####
cl <- getCommandLineArgs()

barcode <- cl$args[1]
ofname <- cl$args[2]

opt <- cl$options
verbose <- opt$verbose
useAnnotMetadata <- !opt$excelMetadata
if(is.null(opt$local)){
  useSynapse <- TRUE
} else {
  useSynapse <- FALSE
  path <- opt$local
}
rawDataVersion <- opt$rawDataVersion

if(useSynapse) synapseLogin()

scriptStartTime <- Sys.time()

if (verbose) message(paste("Processing plate:", barcode, "\n"))

if (useAnnotMetadata) {
  #Build metadata file name list
  if(useSynapse){
    metadataq <- sprintf("select id from syn8466225 WHERE DataType='Metadata' AND Barcode='%s'",barcode)
    metadataTable <- synTableQuery(metadataq)@values
    metadataFiles <- lapply(metadataTable$id, synGet)
    metadataFiles <- list(annotMetadata=getFileLocation(metadataFiles[[1]]))
  } else {
    metadataFiles <- list(annotMetadata=paste0(path,"/",barcode,"_an2omero.csv"))
  }
  
} else {
  metadataFiles <- list(logMetadata = dir(path,pattern = "xml",full.names = TRUE),
                        spotMetadata = dir(path,pattern = "gal",full.names = TRUE),
                        wellMetadata =  dir(path,pattern = "xlsx",full.names = TRUE)
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
  cellDataFilePaths <- dir(paste0(path,"/",rawDataVersion), full.names = TRUE)
  if(length(cellDataFilePaths)==0) stop("No raw data files found")
  dataBWInfo <- data.table(Path=cellDataFilePaths,
                           Well=gsub("_","",str_extract(dir(paste0(path,"/",rawDataVersion)),"_.*_")),
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
                                     DataType = "Quantitative Imaging",
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
