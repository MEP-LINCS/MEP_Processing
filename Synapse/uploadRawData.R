library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(rGithubClient)

synapseLogin()

reportType <- "Preprocessing"
#dataType <- "Quantitative Imaging"
#dataType <- "Metadata"
activityName <- "Annotate data"
#activityName <- "Merge Data and Metadata"
#reportDesc <- "Analysis"


repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="master")
scriptLink <- getPermlink(repo, paste0("MEP-LINCS_",reportType,".R"))

synapseRawDataDir <- "syn5706233"
synapseAnnotatedDataDir <- "syn5706203"
synapseReportDir <- "syn5007815"
synapseMetadataDir <- "syn4997970"

# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  as.list(x %>% select(-filename))
}

# Take row of data frame with filename and annots
# Upload to Synapse and set annotations
uploadToSynapse <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName, 
                  forceVersion=FALSE,
                  executed=thisScript)
  obj
}

#Select datasets to upload
ssDatasets <- rbind(
  data.frame(CellLine=rep(c("PC3"), 4),
             StainingSet=c("SS1", "SS2","SS3","SS2noH3"),
             Preprocess="av1.4",
             Segmentation=c("v2","v2.1","v2.1", "v1"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF7"), 3),
             StainingSet=c("SS1", "SS2","SS3"),
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("YAPC"), 3),
             StainingSet=c("SS1","SS2","SS3"),
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF10A"), 2),
             StainingSet=c("SS1","SS3"),
             Preprocess="av1.4",
             Segmentation=c("v2","v2"),
             stringsAsFactors=FALSE))


getPaths <- function(x){
  dataDir <- paste("/Users/dane/Documents/MEP-LINCS",x$CellLine,x$StainingSet,"RawData/v2", sep = "/")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFiles <- data.frame(filename=list.files(path=dataDir), stringsAsFactors = FALSE) %>%
    mutate(Level=0,
           Consortia="MEP-LINCS",
           DataType="Quantitative Imaging",
           CellLine=x$CellLine,
           StainingSet=x$StainingSet,
           Segmentation=x$Segmentation,
           basename=str_replace(filename, ".*/", "")) %>% 
    mutate(basename=str_replace(basename, "\\.csv", "")) %>% 
    separate(basename, c("Barcode", "Well", "Location"))

  res <- dlply(dataFiles, .(filename), uploadToSynapse, parentId=synapseRawDataDir)
  
}

dataFiles <- do.call(rbind, dlply(ssDatasets[2:nrow(ssDatasets),], c("CellLine","StainingSet"), getPaths))


