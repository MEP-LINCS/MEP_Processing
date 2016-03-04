library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(rGithubClient)

synapseLogin()

repo <- getRepo("MEP-LINCS/MEP_LINCS_Pilot", ref="branch", refName="master")
thisScript <- getPermlink(repo, "uploadRawData.R")

synapseRawDataDir <- "syn5706233"
synapseAnnotatedDataDir <- "syn5706203"

# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  as.list(x %>% select(-filename))
}

# Take row of data frame with filename and annots
# Upload to Synapse and set annotations
uploadToSynapse <- function(x, parentId) {
  browser()
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName="Upload", 
                  forceVersion=False,
                  executed=thisScript)
  obj
}

#Select datasets to upload
ssDatasets <- rbind(
  data.frame(CellLine=rep(c("PC3"), 4),
             StainingSet=c("SS1", "SS2","SS3","SS2noH3"),
             Segmentation=c("v2","v2.1","v2.1", "v1"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF7"), 3),
             StainingSet=c("SS1", "SS2","SS3"),
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("YAPC"), 2),
             StainingSet=c("SS1","SS3"),
             Segmentation=c("v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF10A"), 2),
             StainingSet=c("SS1","SS3"),
             Segmentation=c("v2","v2"),
             stringsAsFactors=FALSE)
)
ssDatasets$fileType <- "tsv"


getPaths <- function(x){
  dataDir <- paste("/Users/dane/Documents/MEP-LINCS",x$CellLine,x$StainingSet,"RawData/v2", sep = "/")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFiles <- data.frame(filename=list.files(path=dataDir, full.names = TRUE), stringsAsFactors = FALSE) %>%
    mutate(level=0,
           CellLine=x$CellLine,
           StainingSet=x$StainingSet,
           Filename=as.character(filename),
           basename=str_replace(filename, ".*/", "")) %>% 
    mutate(basename=str_replace(basename, "\\.csv", "")) %>% 
    separate(basename, c("Barcode", "Well", "Location"))

  res <- dlply(dataFiles[1, ], .(filename), uploadToSynapse, parentId=synapseRawDataDir)
  
}

dataFiles <- do.call(rbind, dlply(ssDatasets, c("CellLine","StainingSet"), getPaths))


