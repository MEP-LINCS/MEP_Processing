library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(rGithubClient)

synapseLogin()

repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="master")
thisScript <- getPermlink(repo, "MEP-LINCS_Preprocessing.R")

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
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName="Normalize and Summarize to Spot Level",
                  contentType="text/tab-separated-values",
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
             stringsAsFactors=FALSE)
)


getPaths <- function(x){
  AllFilePaths <- dir(paste(x$CellLine,x$StainingSet,"AnnotatedData",sep="/"),full.names = TRUE)
  filePaths <- grep(paste0(x$Segmentation,"_",x$Preprocess),AllFilePaths, value=TRUE)
  files <- data.frame(x,filename=filePaths,stringsAsFactors=FALSE, row.names=NULL)
  files$Level <- as.integer(gsub(".*Level|.txt","",files$filename))
  files$Consortia <- "MEP-LINCS"
  files$DataType <- "Quantitative Imaging"

  return(files)
}

dataFiles <- do.call(rbind,dlply(ssDatasets[9,], c("CellLine","StainingSet"), getPaths))

res <- dlply(dataFiles[dataFiles$Level==3&!is.na(dataFiles$Level),], .(filename), uploadToSynapse, parentId=synapseAnnotatedDataDir)
