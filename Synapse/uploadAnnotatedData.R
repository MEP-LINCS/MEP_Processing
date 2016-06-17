library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)

synapseLogin()

scriptName <- "PreprocessingL3L4"
dataType <- "Quantitative Imaging"
#dataType <- "Metadata"
activityName <- "Annotate data"
#activityName <- "Normalize Cell Data"
#reportDesc <- "Analysis"

repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="master")
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",scriptName,".R"))

synapseRawDataDir <- "syn6167751"
synapseAnnotatedDataDir <- "syn5713302"
synapseReportDir <- "syn4939350"
synapseMetadataDir <- "syn5007108"

# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  #as.list(x %>% select(-c(filename,Level1SynID)))
  as.list(x %>% select(-c(filename)))
  
}

# Take row of data frame with filename and annots
# Upload to Synapse and set annotations
uploadToSynapse <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  #used=c(x$Level1SynID),
                  forceVersion=FALSE,
                  executed=scriptLink)
  obj
}


#Select datasets to upload
ssDatasets <- rbind(
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
  data.frame(CellLine=rep(c("MCF10A"), 3),
             StainingSet=c("SS1","SS2","SS3"),
             Drug="none",
             Preprocess="av1.6",
             Segmentation="v2",
             stringsAsFactors=FALSE))

getPaths <- function(x){
  AllFilePaths <- dir("AnnotatedData",full.names = TRUE)
  filePaths <- grep(paste0(x$CellLine,"_",x$StainingSet,"_",x$Drug,"_",x$Segmentation,"_",x$Preprocess),AllFilePaths, value=TRUE)
  files <- data.frame(x,filename=filePaths,stringsAsFactors=FALSE, row.names=NULL)
  files$Level <- as.integer(gsub(".*Level|.txt","",files$filename))
  files$Consortia <- "MEP-LINCS"
  files$DataType <- "Quantitative Imaging"
  #files$Level1SynID<-x[["1"]]
  return(files)
}

dataFiles <- do.call(rbind,dlply(ssDatasets[11:13,], c("CellLine","StainingSet"), getPaths))

res <- dlply(dataFiles[dataFiles$Level==3&!is.na(dataFiles$Level),], .(filename), uploadToSynapse, parentId=synapseAnnotatedDataDir)
