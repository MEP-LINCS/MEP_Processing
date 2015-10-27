library(synapseClient)
library(tidyr)
library(stringr)
library(rGithubClient)

synapseLogin()

repo <- getRepo("kdaily/MEP_LINCS_Pilot", ref="branch", refName="usesynapse")
thisScript <- getPermlink(repo, "uploadRawData.R")

dataDir <- "/home/kdaily/Projects/LINCS/data/Pilot/Raw_Data"
synapseRawDataDir <- "syn4624343"

# Take file names and turn into basic annotation set
# Replace this with a better way to get basic annotations from 
# a standardized source
dataFiles <- data.frame(filename=list.files(path=dataDir, full.names = TRUE)) %>%
  mutate(level=0,
         filename=as.character(filename),
         basename=str_replace(filename, ".*/", "")) %>% 
  mutate(basename=str_replace(basename, "\\.txt", "")) %>% 
  separate(basename, c("Barcode", "Well", "number", "CellLine", "ls", "StainingSet", "index", "dataSubType")) %>% 
  separate(Well, c("row", "column"), sep=1, remove=FALSE) %>%
  mutate(StainingSet=toupper(StainingSet))

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
  obj <- synStore(obj, forceVersion=FALSE, 
                  activityName="Upload", executed=thisScript)
  obj
}

res <- dlply(dataFiles[1:2, ], .(filename), uploadToSynapse, parentId=synapseRawDataDir)
