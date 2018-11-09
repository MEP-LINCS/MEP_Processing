#!/usr/bin/env Rscript

#author: "Mark Dane"

library(MEMA)
library(parallel)
library(stringr)
library(tidyverse)
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(optparse))
library(githubr)


# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-i", "--inputPath"), type="character", default=NULL, metavar="PATH",
                help="Path to local input data directory or Synapse ID for a File View."),
    make_option(c("--synapseStore"), type="character", default=NULL, metavar="SYNAPSEID",
                help="Store output file in Synapse directory (provide Synapse ID of Folder to store).")
  )
  parser <- OptionParser(usage = "%prog [options] BARCODE OUTPUTFILE", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#cl <- getCommandLineArgs()
barcode <- cl$args[1]
ofname <- cl$args[2]
opt <- cl$options
verbose <- opt$verbose
if(file.exists(opt$inputPath)){
  useSynapse <- FALSE
} else {
  useSynapse <- TRUE
}

if (verbose) message(sprintf("Summarizing cell to spot data for plate %s", barcode))
functionStartTime<- Sys.time()
startTime<- Sys.time()
seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")

#Read in the plate's cell level data
if (useSynapse) {
  suppressMessages(synapseLogin())
  level <- "1"
  levelQuery <- sprintf('SELECT id,Segmentation,Preprocess,DataType,Study,Consortia,StainingSet,CellLine,Drug from %s WHERE Barcode="%s" AND Level="%s"',
                        opt$inputPath, barcode, level)
  levelRes <- synTableQuery(levelQuery)
  
  if (nrow(levelRes@values) > 1) {
    stop(sprintf("Found more than one Level 1 file for barcode %s", barcode))
  }
  
  dataPath <- getFileLocation(synGet(levelRes@values$id))
  
  imageIdQuery <- sprintf("SELECT id from %s WHERE Barcode='%s' AND DataType='ImageID'",
                          opt$inputPath, barcode)
  imageIdRes <- synTableQuery(imageIdQuery)
  
  clarionIdQuery <- sprintf('SELECT id from %s WHERE Barcode="%s" AND DataType="ClarionID"',
                            opt$inputPath, barcode)
  clarionIdRes <- synTableQuery(clarionIdQuery)
  
  if (nrow(imageIdRes@values) == 1) imageIdPath <- getFileLocation(synGet(imageIdRes@values$id))

  if (nrow(clarionIdRes@values) == 1)   clarionIdPath <- getFileLocation(synGet(clarionIdRes@values$id))
  
  
} else {
  dataPath <- paste0(opt$inputPath, "/",barcode, "_Level1.tsv")
  imageIdPath <- paste0(opt$inputPath, "/",barcode, "_imageIDs.tsv")
  clarionIdPath <- paste0(opt$inputPath, "/",barcode, "_clarionIDs.tsv")
}
######Debug   move this into MEMA package
#' Read in and merge the Omero URLs
#' 
#' Adds Omero image IDs based on the WellName values
#' 
#' @param path The path to the file named barcode_imageIDs.tsv
#' @return a datatable with WellIndex, ArrayRow, ArrayColumn and ImageID columns
#' @export
getOmeroIDs <- function(path){
  dt <- fread(path)[,list(WellName,Row,Column,ImageID)]
  # if(any(grepl("LI9",unique(dt$WellName)))){
  #   #Extract well index and convert to an integer in a 96 well plate
  #   well <- gsub("_.*","",gsub(".*_Well","",dt$WellName))
  #   wellRow <- gsub("[[:digit:]]*","",well)
  #   wellColumn <- as.integer(gsub("[[:alpha:]]","",well))
  #   dt$WellIndex <- (match(wellRow,LETTERS)-1)*12+wellColumn
  # } else {
    #Extract well index and convert to an integer in an 8 well plate
    dt <- dt[,WellIndex := as.integer(gsub(".*_Well","",WellName))]
  # }

  setnames(dt,"Row","ArrayRow")
  setnames(dt,"Column","ArrayColumn")
  dt[,WellName := NULL]
  return(dt)
}
#####
cDT <- fread(dataPath)
if(exists("imageIdPath")){
if(file.exists(imageIdPath)) omeroIDs <- getOmeroIDs(imageIdPath)
}



if(exists("clarionIdPath")){
  if(file.exists(clarionIdPath)){
    clarionIDs <- getOmeroIDs(clarionIdPath)
    setnames(clarionIDs, "ImageID", "ClarionID")
    if (useSynapse) levelRes@values$Preprocess <- paste0(levelRes@values$Preprocess,".1")
  } 
}

#Count the cells at each spot at the cell level as needed by createl3
cDT <- cDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]

#Add the DAPI intensities in all nuclei
cDT <- cDT[,Spot_PA_TotalDapiIntensity := sum(Nuclei_CP_Intensity_IntegratedIntensity_Dapi),by="Barcode,Well,Spot"]

#Add the EdU intensities in all nuclei if present
if("Nuclei_CP_Intensity_IntegratedIntensity_EdU" %in% colnames(cDT)) cDT <- cDT[,Spot_PA_TotalEdUIntensity := sum(Nuclei_CP_Intensity_IntegratedIntensity_EdU),by="Barcode,Well,Spot"]

#Add a luminal count feature
if("Cytoplasm_PA_Gated_KRT19Positive" %in% colnames(cDT)) cDT <- cDT[,Cytoplasm_PA_Gated_KRT19Positive_SpotCellCount := sum(Cytoplasm_PA_Gated_KRT19Positive),by="Barcode,Well,Spot"]

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

# create an organization feature
calcOrganization <- function(df){
#Use the existing outer cell labels which divide the dataset into 50% interior
  df %>%
    dplyr::group_by(Well, Spot) %>%
    dplyr::mutate(Spot_PA_InteriorKRT19HighProportion = sum(Cytoplasm_PA_Gated_KRT19Positive&!Spot_PA_OuterCell)/sum(!Spot_PA_OuterCell),
           Spot_PA_PerimeterKRT19LowProportion = sum(!Cytoplasm_PA_Gated_KRT19Positive&Spot_PA_OuterCell)/sum(Spot_PA_OuterCell),
           Spot_PA_KRT19Organization = (Spot_PA_InteriorKRT19HighProportion+Spot_PA_PerimeterKRT19LowProportion)/2) %>%
    data.table()
}

#cDT <- calcOrganization(cDT)


#median summarize the rest of the signals to the spot level
signals <- summarizeToSpot(cDT, seNames)

if(exists("proportions")) {
  spotDT <- merge(signals,proportions)
} else {
  spotDT <- signals
}

#Set a threshold for the loess well level QA Scores
if(sum(c("ArrayRow","ArrayColumn") %in% colnames(spotDT))==2) spotDT <- QASpotData(spotDT, lthresh = .6)

#Merge in Omero imageID links
if(exists("omeroIDs")) spotDT <- merge(spotDT, omeroIDs,
                                       by=c("WellIndex","ArrayRow","ArrayColumn"))
#Merge in Clarion imageID links
if(exists("clarionIDs")) spotDT <- merge(spotDT, clarionIDs,
                                         by=c("WellIndex","ArrayRow","ArrayColumn"))

#Develop load roboust data if it exists
#If there is robust data, load, clean and create compatable well values
if(!is_empty(dir(paste0(cl[["options"]]$inputPath,"/QA"), full.names = TRUE))){
  
  rb <- map(dir(paste0(cl[["options"]]$inputPath,"/QA"), full.names = TRUE), read_csv) %>% bind_rows() %>%
    mutate(Spot = imageID,
           Well = str_extract(imageName, "Well."),
           Well = str_remove(Well, "Well"),
           Well = as.integer(Well),
           Well = c("A01","A02","A03","A04","B01","B02","B03","B04")[Well]) %>%
    select(-X25) %>%
    purrr::set_names(~str_replace_all(.," ", "_"))
  
  #Merge roboust and CP data
  spotDT <- spotDT %>%
    full_join(rb, by=c("Well","Spot"))
} 
  
if(verbose) message("Writing spot level data")
writeTime<-Sys.time()
#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}
for (j in colnames(spotDT)) data.table::set(spotDT, j = j, value = shorten(spotDT[[j]]))

fwrite(spotDT, file=ofname, sep = "\t", quote=FALSE)

if(!is.null(cl$options$synapseStore)){
  if(verbose) message(sprintf("Writing to Synapse Folder %s", opt$synapseStore))
  #get permlink from GitHub
  scriptLink <- "https://github.com/MEP-LINCS/MEP_Processing/"
  repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)
  if(!class(repo)=="try-error" ) scriptLink <- getPermlink(repo, "Pipeline/PreprocessMEMALevel2.R")
  synFile <- File(ofname, parentId=opt$synapseStore)
  synSetAnnotations(synFile) <- list(CellLine = levelRes@values$CellLine,
                                     Barcode = barcode,
                                     Study = levelRes@values$Study,
                                     Preprocess = levelRes@values$Preprocess,
                                     DataType = levelRes@values$DataType,
                                     Consortia = levelRes@values$Consortia,
                                     Drug = levelRes@values$Drug,
                                     Segmentation = levelRes@values$Segmentation,
                                     StainingSet = levelRes@values$StainingSet,
                                     Level = "2")
  
  synFile <- synStore(synFile,
                      used=c(levelRes@values$id, imageIdRes@values$id),
                      executed=scriptLink,
                      forceVersion=FALSE)
}

if(verbose) message(paste("Write time:", Sys.time()-writeTime,"\n"))
if(verbose) message(paste("Elapsed time:", Sys.time()-functionStartTime, "\n"))

write.csv(colnames(spotDT),file = "../QAFeatureNames")

