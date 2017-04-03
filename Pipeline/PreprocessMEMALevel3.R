#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2/2017

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(RUVnormalize)
library(ruv)
library(stringr)
library(tidyr)
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-r", "--rawDataVersion"), type="character", default=NULL,
                help="Raw data version [default \"%default\"]"),
    make_option(c("-k", "--k"), type="integer", default=256,
                help="Number of factors to use in RUV normalization [default %default]",
                metavar="number"),
    make_option(c("-i", "--inputPath"), type="character", default=NULL, metavar="PATH",
                help="Path to local input data directory or Synapse ID for a File View."),
    make_option(c("--synapseStore"), type="character", default=NULL, metavar="SYNAPSEID",
                help="Store output file in Synapse directory (provide Synapse ID of Folder to store).")
  )
  parser <- OptionParser(usage = "%prog [options] STUDY OUTPUTFILE", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
###Debug
cl <-list(options=list(verbose=TRUE,
                       local="/lincs/share/lincs_user",
                       k=256),
          args=c("HMEC122L_SS1",
                 "/lincs/share/lincs_user/study/HMEC122L_SS1/Annotated/HMEC122L_SS1_Level3.tsv")
)
####
cl <- getCommandLineArgs()

studyName <- cl$args[1]
ofname <- cl$args[2]

opt <- cl$options
verbose <- opt$verbose

if(file.exists(opt$inputPath)){
  useSynapse <- FALSE
} else {
  useSynapse <- TRUE
}

k <- opt$k

startTime <- Sys.time()

#Read the annotated data for all plates in the study
if(useSynapse){
  suppressPackageStartupMessages(library(synapseClient))
  suppressMessages(synapseLogin())
  level <- "2"
  levelQuery <- sprintf('SELECT id,rawDataVersion from %s WHERE Study="%s" AND Level="%s"',
                        opt$inputPath, studyName, level)
  levelRes <- synTableQuery(levelQuery)
  
  if (nrow(levelRes@values) > 1) {
    stop(sprintf("Found more than one Level %s file for Study %s: %s", 
                 level, studyName, paste(levelRes@values$id, collapse=",")))
  }
  
  rawDataVersion <- levelRes@values$rawDataVersion
  
  dataPath <- getFileLocation(synGet(levelRes@values$id))
  
} else {
  dataPath <- opt$inputPath
  rawDataVersion <- opt$rawDataVersion
}

slDT <- getSpotLevelData(studyName, dataPath)

signalsMinimalMetadata <- grep("_SE",
                               grep("_CP_|_PA_|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^Drug1Conc$|^ArrayRow$|^ArrayColumn$|^CellLine$",
                                    colnames(slDT), value=TRUE), 
                               value=TRUE, invert=TRUE)

#RUVLoess normalize all signals
if(!k==0){
  if(verbose)  message(paste("Normalizing", studyName,"\n"))
  #slDT <- slDT[!is.na(slDT$Cytoplasm_CP_AreaShape_Compactness)&!is.na(slDT$Cytoplasm_CP_AreaShape_Eccentricity),]
  nDT <- normRUVLoessResiduals(slDT[,signalsMinimalMetadata, with = FALSE], k)
  nDT$NormMethod <- "RUVLoessResiduals"
  slDT$k <- k
  slDT <- merge(slDT, nDT, by = c("BW","PrintSpot"))
} else {
  slDT$NormMethod <- "none"
  slDT$k <- k
}

#Add QA flags to the data
slDT <- QASpotLevelData(slDT, lowSpotCellCountThreshold=5,
                        lowRegionCellCountThreshold = 0.4,
                        lowWellQAThreshold = .7)

if(verbose) message(paste("Writing level 3 file to disk\n"))
fwrite(data.table(slDT), file=ofname, sep = "\t", quote=FALSE)

if(!is.null(opt$synapseStore)) {
  if(verbose) message(sprintf("Writing to Synapse Folder %s", opt$synapseStore))
  synFile <- File(ofname, parentId=opt$synapseStore)
  synSetAnnotations(synFile) <- list(CellLine = unique(slDT$CellLine),
                                     # Barcode = barcode,
                                     Study = unique(slDT$Study),
                                     Preprocess = "v1.8",
                                     DataType = "Quanititative Imaging",
                                     Consortia = "MEP-LINCS",
                                     Drug = unique(slDT$Drug),
                                     Segmentation = rawDataVersion,
                                     StainingSet = gsub("Layout.*","",unique(slDT$StainingSet)),
                                     Level = "3")
  
  synFile <- synStore(synFile,
                      used=c(levelRes@values$id),
                      forceVersion=FALSE)
}

message(paste("Elapsed time to normalize ",studyName, Sys.time()-startTime, "\n"))


