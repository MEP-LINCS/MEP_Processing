#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
library(tidyverse)
library(synapseClient)
library(MEMA)
library(parallel)
library(stringr)
library(readr)
suppressPackageStartupMessages(library(dplyr))
library(optparse)
library(githubr)


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

# #Specify the command line options
# cl <- getCommandLineArgs()
# barcode <- cl$args[1]
# ofname <- cl$args[2]
# opt <- cl$options
# verbose <- opt$verbose
# useAnnotMetadata <- !opt$excelMetadata

#Specify the command line options
if(!interactive()){
  cl <- getCommandLineArgs()
  barcode <- cl$args[1]
  ofname <- cl$args[2]
  opt <- cl$options
  verbose <- opt$verbose
  useAnnotMetadata <- !opt$excelMetadata
  inputPath <- opt$inputPath
  rawDataVersion <- opt$rawDataVersion
  synpaseStore <- opt$synapseStore
} else {
  barcode <- "LI8V01180a"
  ofname <- "/lincs/share/lincs_user/LI8V01180a/Analysis/LI8V01180a_Level1.tsv"
  verbose <- FALSE
  dataDir <- "/lincs/share/lincs_user/incellSlides"
  useAnnotMetadata <- TRUE
  inputPath <- "/lincs/share/lincs_user/LI8V01180a/Analysis"
  rawDataVersion <- "v2"
  synapseStore <- NULL
}

  
if(file.exists(inputPath)){
  useSynapse <- FALSE
} else {
  useSynapse <- TRUE
}

rawDataVersion <- opt$rawDataVersion

if(useSynapse) synapseLogin()

scriptStartTime <- Sys.time()

if (verbose) message(paste("Processing plate:", barcode, "\n"))

if (useAnnotMetadata) {
  #Build metadata file name list
  if(useSynapse){
    metadataq <- sprintf("select id from %s WHERE DataType='Metadata' AND Barcode='%s'",
                         inputPath, barcode)
    metadataTable <- synTableQuery(metadataq)@values
    metadataFiles <- lapply(metadataTable$id, synGet)
    metadataFiles <- list(annotMetadata=getFileLocation(metadataFiles[[1]]))
  } else {
    metadataFiles <- list(annotMetadata=paste0(inputPath,"/",barcode,"_an2omero.csv"))
  }
  
} else {
  metadataFiles <- list(logMetadata = dir(inputPath,
                                          pattern = "xml",full.names = TRUE),
                        spotMetadata = dir(inputPath,
                                           pattern = "gal",full.names = TRUE),
                        wellMetadata =  dir(inputPath,
                                            pattern = "xlsx",full.names = TRUE)
  )
}

#Get all metadata
metadata <- getMetadata(metadataFiles, useAnnotMetadata)

#reimaged plates may have a sufgix on the data that is not inlcuded in the 
#barcode value in the metadata. Add the suffix to the metadata barcode
suffix <- str_extract(barcode, "[[:alpha:]]$")
if(!is_empty(suffix)) {
  metadata$Barcode <- paste0(metadata$Barcode,suffix)
}

#Get image quality data if available
QA <- lapply(list.files(paste0(inputPath,"/v2"),pattern="SummaryImageData",full.names = TRUE), function(x) suppressWarnings(read_csv(file=x, col_types = cols()))) %>%
  bind_rows %>%
  select(-matches("^CP |^X"))
if(!nrow(QA)==0){
  colnames(QA) <- str_replace_all(colnames(QA)," ","_")
  QA <- rename(QA,Spot=imageID, Well=imageGroupName)
  QA$Well[seq(from=700,to=5600,by=700)] <- QA$Well[seq(from=699,to=5600,by=700)]
  
  QA <- mutate(QA, Barcode = str_extract(imageName,".*?_"))
  QA$Barcode <- str_replace_all(QA$Barcode, "_|Plate","")
}

#Gather filenames and metadata of level 0 files
if(useSynapse){
  q <- sprintf("select id,Barcode,Level,Well,StainingSet,Location,Study from %s WHERE Level='0' AND Barcode='%s'",
               inputPath, barcode)
  rawFiles <- synTableQuery(q)
  dataBWInfo <- rawFiles@values
  
  # Download raw files, or get from cache if already downloaded
  res <- lapply(dataBWInfo$id, synGet)
  
  # Get the path on disk to each file
  cellDataFilePaths <- unlist(lapply(res, getFileLocation))
  dataBWInfo$Path <- cellDataFilePaths
} else {
  cellDataFilePaths <- dir(paste0(inputPath,"/",rawDataVersion), pattern="Nuclei|Cells|Cytoplasm",full.names = TRUE)
  if(length(cellDataFilePaths)==0) stop("No raw data files found")
  dataBWInfo <- data.table(Path=cellDataFilePaths,
                           Well=gsub("_","",str_extract(dir(paste0(inputPath,"/",rawDataVersion),pattern="Nuclei|Cells|Cytoplasm"),"_.*_")),
                           Location=str_extract(cellDataFilePaths,"Nuclei|Cytoplasm|Cells|Image"))
}
if(length(cellDataFilePaths) == 0) stop("No raw data files found")

#Gather data from either Cell Profiler or INCell
if("Nuclei" %in% dataBWInfo$Location) { #Cell Profiler data
  dtL <- getCPData(dataBWInfo = dataBWInfo, verbose=verbose, curatedOnly = FALSE)
} else if (any(grepl("96well", dataBWInfo$Path))) { #GE INCell data
  dtL <- getICData(cellDataFilePaths = cellDataFilePaths, 
                   endPoint488 = unique(metadata$Endpoint488),
                   endPoint555 = unique(metadata$Endpoint555),
                   endPoint647 = unique(metadata$Endpoint647),
                   verbose=verbose)
} else {
  stop("Only data from Cell Profiler or INCell pipelines are supported")
}

#Detect 96 well CP data and assign well and spot values
if(any(grepl("Row",dtL[[1]]$Well,ignore.case = TRUE))){
  dtL <- lapply(dtL, formatCP96Well,
                nrArrayRows= max(metadata$ArrayRow),
                nrArrayColumns = max(metadata$ArrayColumn))
} 

#Clean up legacy issues in column names and some values
dtL <- lapply(dtL, cleanLegacyIssues)

# Filter our debris and cell clusters
dtL <- lapply(dtL, function(dt){
  filterObjects(dt,nuclearAreaThresh = 200, nuclearAreaHiThresh = 4000)})

# Add local XY and polar coordinates
# dtL <- mclapply(dtL, addPolarCoords, mc.cores = detectCores())
dtL <- lapply(dtL, addPolarCoords)

# Add spot level normalizations for median intensities
# dtL <- mclapply(dtL,spotNormIntensities, mc.cores = detectCores())
dtL <- lapply(dtL,spotNormIntensities)

# Add adjacency values
# dtL <- mclapply(dtL, calcAdjacency, mc.cores = detectCores())
dtL <- lapply(dtL, calcAdjacency)

if(nrow(QA)>0){
  # Merge the data with metadata and image QA results
  cDT <- merge(rbindlist(dtL),metadata,by=c("Barcode","Well","Spot"),all = TRUE) %>%
    merge(QA,by=c("Barcode","Well","Spot"),all = TRUE)
} else {
  # Merge the data with metadata and image QA results
  cDT <- merge(rbindlist(dtL),metadata,by=c("Barcode","Well","Spot"),all = TRUE)
}

#####
#Temp filter for 0 area donuts
if(exists("cDT$Cytoplasm_CP_AreaShape_Area"))
cDT <- cDT[!is.na(Cytoplasm_CP_AreaShape_Area),]
if("Nuclei_CP_AreaShape_Area" %in% colnames(cDT))
  cDT <- cDT[!is.na(Nuclei_CP_AreaShape_Area),]
#####
# Gate cells where possible
cDT <- gateCells(cDT)

#Remove Parent features
cDT <- select(cDT,-contains("Parent"))

#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}
for (j in colnames(cDT)) data.table::set(cDT, j = j, value = shorten(cDT[[j]]))

if(verbose) message("Writing cell level data\n")
fwrite(cDT, file=ofname, sep = "\t", quote = FALSE)
if(!is.null(synapseStore)){
  #get permlink from GitHub
  scriptLink <- "https://github.com/MEP-LINCS/MEP_Processing/"
  repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)
  if(!class(repo)=="try-error" ) scriptLink <- getPermlink(repo, "Pipeline/PreprocessMEMACell.R")
  synFile <- File(ofname, parentId=synapseStore)
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
                      executed=scriptLink,
                      forceVersion=FALSE)
}

if(verbose) message(paste("Elapsed time:", Sys.time()-scriptStartTime, "\n"))

