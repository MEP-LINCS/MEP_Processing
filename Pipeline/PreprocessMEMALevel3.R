#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2/2017

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(RUVnormalize)
library(ruv)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
library(githubr)

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
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

message("Preprocessing level 3 data ")

#Specify the command line options
if(!interactive()){
  cl <- getCommandLineArgs()
  studyName <- cl$args[1]
  ofname <- cl$args[2]
  inputPath <- cl$options$inputPath
  verbose <- cl$options$verbose
  k <- cl$options$k
  synapseStore <- cl$options$synapseStore
} else {
  studyName <- "cama1_highserum_fulvest_mema"
  ofname <- "/lincs/share/lincs_user/study/cama1_highserum_fulvest_mema/Annotated/cama1_highserum_fulvest_mema_Level3.tsv"
  inputPath <- "/lincs/share/lincs_user"
  verbose <- FALSE
  inputPath <- "/lincs/share/lincs_user"
  k <- 256
  synapseStore <- FALSE
}

message(studyName)

if(file.exists(inputPath)){
  useSynapse <- FALSE
} else {
  useSynapse <- TRUE
}

startTime <- Sys.time()

#Read the annotated data for all plates in the study
if(useSynapse){
  suppressPackageStartupMessages(library(synapseClient))
  suppressMessages(synapseLogin())
  level <- "2"
  levelQuery <- sprintf('SELECT id,Segmentation,Preprocess,Study,DataType,Consortia,StainingSet,CellLine,Drug from %s WHERE Study="%s" AND Level="%s"',
                        inputPath, studyName, level)
  levelRes <- synTableQuery(levelQuery)
  
  dataPaths <- lapply(levelRes@values$id,
                      function(x) getFileLocation(synGet(x)))
  
} else {
  barcodes <- getBarcodes(studyName, synId = "syn10846457")
  dataPaths <- barcodes %>%
    lapply(function(barcode, path){
      paste0(path,"/",barcode,"/Analysis/",barcode,"_Level2.tsv")
    }, path=inputPath)
}

#Special handling for MEP analysis, remove NID1 and ELN data
slDT <- getSpotLevelData(dataPaths) %>%
  # filter(!ECMp=="NID1|1",
  #        !ECMp=="ELN|3",
  #        !Ligand %in% c("KNG1|HMW", "LYVE1", "THPO|1", "JAG2|Long"),
  #        !(Barcode=="LI8X00652"&Ligand=="FBS"),
  #        !(Barcode=="LI8X00656"&Ligand=="FBS")) %>%
  data.table()


slDT$Drug[slDT$Drug== "air"] <- "DMSO"
slDT$Drug1[slDT$Drug1== "air_Own"] <- "DMSO_chebi28262"

signalsMinimalMetadata <- grep("_SE",
                               grep("_CP_|_PA_|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^Drug1Conc$|^ArrayRow$|^ArrayColumn$|^CellLine$",
                                    colnames(slDT), value=TRUE), 
                               value=TRUE, invert=TRUE)

#RUVLoess normalize all signals
if(!k==0){
  if(verbose)  message(paste("Normalizing", studyName,"\n"))
  slDT <- slDT[!is.na(slDT$Cytoplasm_CP_AreaShape_Compactness)&!is.na(slDT$Cytoplasm_CP_AreaShape_Eccentricity),]
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

#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}
for (j in colnames(slDT)) data.table::set(slDT, j = j, value = shorten(slDT[[j]]))

if(verbose) message(paste("Writing level 3 file to disk\n"))
fwrite(data.table(slDT), file=ofname, sep = "\t", quote=FALSE)

if(!is.null(synapseStore)) {
  if(verbose) message(sprintf("Writing to Synapse Folder %s", synapseStore))
  #get permlink from GitHub
  scriptLink <- "https://github.com/MEP-LINCS/MEP_Processing/"
  repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)
  if(!class(repo)=="try-error" ) scriptLink <- getPermlink(repo, "Pipeline/PreprocessMEMALevel3.R")
  synFile <- File(ofname, parentId=synapseStore)
  synSetAnnotations(synFile) <- list(CellLine = unique(slDT$CellLine),
                                     Study = unique(levelRes@values$Study),
                                     Preprocess = unique(levelRes@values$Preprocess),
                                     DataType = unique(levelRes@values$DataType),
                                     Consortia = unique(levelRes@values$Consortia),
                                     Drug = unique(levelRes@values$Drug),
                                     Segmentation = unique(levelRes@values$Segmentation),
                                     StainingSet = unique(levelRes@values$StainingSet),
                                     Level = "3")
  
  synFile <- synStore(synFile,
                      used=c(levelRes@values$id),
                      executed=scriptLink,
                      forceVersion=FALSE)
}

message(paste("Elapsed time to normalize ",studyName, Sys.time()-startTime, "\n"))


