library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)
library(MEMA)

synapseLogin()

reportType <- "Preprocessing"
dataType <- "ImageIDs"
activityName <- "Annotate data"

repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="master")

synapseMetadataDir <- "syn7213947"

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
                  forceVersion=FALSE)
  obj
}


createAnnotations <- function(studyName, path) {
  #Get barcodes in the study
  barcodes <- getBarcodes(studyName)
  
  #create a dataframe with path, filename, and annotations
  metadataFiles <- data.frame(filename = paste0(path,"/",barcodes,"/Analysis/",barcodes,"_imageIDs.tsv"),
                              Consortia="MEP-LINCS",
                              DataType="ImageID",
                              CellLine=gsub("_.*","",studyName),
                              StainingSet=gsub(".*_","",studyName),
                              Barcode=barcodes,
                              Study = studyName,
                              stringsAsFactors = FALSE
  )
}

metadataFiles <- createAnnotations(studyName = "MCF10A_SS3", path = "/lincs/share/lincs_user")
res <- dlply(metadataFiles, .(filename), uploadToSynapse, parentId=synapseMetadataDir)

for(study in c("HMEC240L_SS1","HMEC240L_SS4","HMEC122L_SS1","HMEC122L_SS4")){
  metadataFiles <- createAnnotations(studyName = study, path = "/lincs/share/lincs_user")
  res <- dlply(metadataFiles, .(filename), uploadToSynapse, parentId=synapseMetadataDir)
}

