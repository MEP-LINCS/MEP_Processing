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
  dataFile <- data.frame(DatasetName=x$DatasetName,
                         CellLine=x$CellLine,
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
synapseRawDataDir <- "syn7804139" #MEMA Experiments project

#Select datasets to upload
datasetAnns <- rbind(
  data.frame(DatasetName="",
             CellLine=rep(c("PC3"), 4),
             StainingSet=c("SS1", "SS2","SS3","SS2noH3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2.1","v2.1", "v1"),
             stringsAsFactors=FALSE),
  data.frame(DatasetName="",
             CellLine=rep(c("MCF7"), 3),
             StainingSet=c("SS1", "SS2","SS3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(DatasetName="",
             CellLine=rep(c("YAPC"), 3),
             StainingSet=c("SS1","SS2","SS3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(DatasetName="",
             CellLine=rep(c("MCF10A"), 4),
             StainingSet=c("SS1","SS2","SS3","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(DatasetName="",
             CellLine=rep(c("HMEC240L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(DatasetName="",
             CellLine=rep(c("HMEC122L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(DatasetName="MCF10AHighRep3",
             CellLine="MCF10A",
             StainingSet="SS4",
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
datasetName=c("MCF10AHighRep3")
#Get all files in all selected datasets
fileNamesList <- apply(datasetBarcodes[datasetBarcodes$DatasetName %in% datasetName,], 1, getMEMADataFileNames)
fileNames <- rbindlist(fileNamesList)
fileNames$DatasetName <- datasetName
#Filter out metadata files
fileNames <- fileNames[grepl("Raw",fileNames$Type),]
#Add in the dataset annotations
ssDatasets <- merge(datasetAnns, fileNames, by=c("DatasetName"))

#Upload the level o files without a provenance to lower level files
res <- dlply(ssDatasets[ssDatasets$Barcode %in% unique(ssDatasets$Barcode)[4],], c("Path"), uploadFileLevel0)
