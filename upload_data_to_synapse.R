library(synapser)
library(tidyverse)
library(purrrlyr)
library(githubr)

synLogin()

reportType <- "Preprocessing"
dataType <- "metadata"
activityName <- "Summarize data"

#repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)

MDD_MEMA_Data <- "syn22264760" 
MDD_MEMA_data_file_suffix <- "_Level2.tsv"
MDD_MEMA_Metadata <- "syn22264761"
MDD_MEMA_metadata_file_suffix <- "_an2omero.csv"
synapse_destination <- MDD_MEMA_Metadata
file_suffix <- MDD_MEMA_metadata_file_suffix

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
  obj_stored <- synStore(obj, 
                         activityName=activityName, 
                         forceVersion=FALSE)
  res <- synSetAnnotations(obj_stored, annots)

  obj_stored
}


createAnnotations <- function(studyName, path) {
  #Get barcodes in the study
  barcodes <- MEMA::getBarcodes(studyName)
  
  #create a dataframe with path, filename, and annotations
  metadataFiles <- tibble(filename = paste0(path,"/",barcodes,"/Analysis/",barcodes,file_suffix),
                          Drug = "none",
                          Level = "",
                          Study = studyName,
                          Barcode=barcodes,
                          CellLine="MCF10A",
                          DataType="Metadata",
                          StainingSet = "SSF",
                          Consortia="MEP-LINCS",
                          segmentation = "",
                          assay = "MEMA"
  )
}

metadataFiles <- createAnnotations(studyName = "mcf10a_egf_ssf", path = "/lincs/share/lincs_user")
res <- metadataFiles %>%
 by_row(uploadToSynapse,  parentId = synapse_destination)

# 
# scriptLink <- "https://github.com/MEP-LINCS/MEP_Processing/"
# repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)
# if(!class(repo)=="try-error" ) scriptLink <- getPermlink(repo, "Pipeline/PreprocessMEMACell.R")
# synFile <- File(ofname, parentId=synapseStore)
# synSetAnnotations(synFile) <- list(CellLine = unique(cDT$CellLine),
#                                    Barcode = barcode,
#                                    Study = unique(cDT$Study),
#                                    Preprocess = "v1.8",
#                                    DataType = "Quantitative Imaging",
#                                    Consortia = "MEP-LINCS",
#                                    Drug = unique(cDT$Drug),
#                                    Segmentation = rawDataVersion,
#                                    StainingSet = gsub("Layout.*","",unique(cDT$StainingSet)),
#                                    Level = "1")
# 
# synFile <- synStore(synFile,
#                     used=c(dataBWInfo$id, metadataTable$id),
#                     executed=scriptLink,
#                     forceVersion=FALSE)
# }

