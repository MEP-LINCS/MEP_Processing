library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)
library(limma)
library(XLConnect)
library(data.table)

source("MEP_LINCS/Release/MEPLINCSFunctions.R")

# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  if("Level3ID" %in% colnames(x)) x <- select(x,-Level3ID)
  if("Used3" %in% colnames(x)) x <- select(x,-Used3)
  if("Used4" %in% colnames(x)) x <- select(x,-Used4)
  if("Level1ID" %in% colnames(x)) x <- select(x,-Level1ID)
  as.list(x %>% select(-c(filename)))
}

uploadToSynapseLevel0 <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  forceVersion=FALSE)
  obj
}

uploadFileLevel0 <- function(x){
  dataFile <- data.frame(CellLine=x$CellLine,
                         StainingSet=x$StainingSet,
                         Drug=x$Drug,
                         Preprocess=x$Preprocess,
                         Segmentation=x$Segmentation,
                         Barcode=x$Barcode,
                         Well=x$Well,
                         Location=x$Location,
                         DataType=dataType,
                         Consortia="MEP-LINCS",
                         Level="0",
                         filename=x$Path,
                         stringsAsFactors = FALSE)
  uploadToSynapseLevel0(dataFile, parentId=synapseRawDataDir)
}


synapseLogin()
synapseRawDataDir <- "syn7525205"

#Select datasets to upload
datasetAnns <- rbind(
  data.frame(CellLine=rep(c("PC3"), 4),
             StainingSet=c("SS1", "SS2","SS3","SS2noH3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2.1","v2.1", "v1"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF7"), 3),
             StainingSet=c("SS1", "SS2","SS3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("YAPC"), 3),
             StainingSet=c("SS1","SS2","SS3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF10A"), 4),
             StainingSet=c("SS1","SS2","SS3","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC240L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC122L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE))

#Begin data or report specific code
#scriptName <- "PreprocessingL1L2"
dataType <- "Quantitative Imaging"
activityName <- "Analyze Images"

# #Get GitHub link for next uploads
#Read the dataset manifest 
datasetBarcodes <- readWorksheetFromFile("DatasetManifest.xlsx", sheet=1)
#Get the file names of all files in these datasets
datasetName=c("MCF10A_SS1","MCF10A_SS2","MCF10A_SS3","HMEC240L_SS1", "HMEC240L_SS4","HMEC122L_SS4","HMEC122L_SS1")
#Get all files in all selected datasets
fileNamesList <- apply(datasetBarcodes[datasetBarcodes$DatasetName %in% datasetName,], 1, getMEMADataFileNames)
#Name the file sets by their dataset names
names(fileNamesList) <-datasetBarcodes$DatasetName[datasetBarcodes$DatasetName %in% datasetName]
#Use the dataset names to add cell line and staining set columns
fileNamesList <- lapply(datasetBarcodes$DatasetName[datasetBarcodes$DatasetName %in% datasetName],function(x){
  ds <- fileNamesList[[x]]
  ds$CellLine <- gsub("_.*","",x)
  ds$StainingSet <- gsub(".*_","",x)
  return(ds)
})
fileNames <- rbindlist(fileNamesList)
#Filter out metadata files
fileNames <- fileNames[grepl("Raw",fileNames$Type),]
#Add in the dataset annotations
ssDatasets <- merge(datasetAnns,fileNames, by=c("CellLine","StainingSet"))

#Upload the level o files without a provenance to lower level files
res <- dlply(ssDatasets, c("Path"), uploadFileLevel0)