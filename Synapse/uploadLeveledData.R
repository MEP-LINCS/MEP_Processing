library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)

#Debug intermediate edits have code failing!!!
synapseLogin()

reportType <- "Preprocessing"
dataType <- "Quantitative Imaging"
#dataType <- "Metadata"
activityName <- "Annotate data"
#activityName <- "Normalize Cell Data"
#reportDesc <- "Analysis"


repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="master")
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",reportType,".R"))

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

#Code to get the level 1 syn ID which may not exist if not uploading
#files <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
#colnames(files)<-gsub("file.","",colnames(files))

#ssDatasets <- merge(ssDatasets,files)
#ssDatasets <- select(ssDatasets, -name)
#ssDatasets <- reshape2::dcast(ssDatasets,CellLine+StainingSet+Drug+Preprocess+Segmentation~Level, value.var="id")

uploadReport <- function(x){
  dataDir <- paste("AnnotatedData/",sep="/")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFile <- data.frame(CellLine=x$CellLine,
                         StainingSet=x$StainingSet,
                         Drug=x$Drug,
                         Preprocess=x$Preprocess,
                         Segmentation=x$Segmentation,
                         DataType=dataType,
                         Consortia="MEP-LINCS",
                         Level=3,
                        # Level1SynID=x[[1]],
                         stringsAsFactors = FALSE)
  
  dataFile$filename <- paste(dataDir,paste0(paste(x$CellLine,x$StainingSet,x$Drug,x$Segmentation,x$Preprocess,"Level1",sep="_"),".txt"),sep="/")
  filename=list.files(path=dataDir), stringsAsFactors = FALSE
  uploadToSynapse(dataFile, parentId=synapseAnnotatedDataDir)
  
}

res <- dlply(ssDatasets[11,], c("CellLine","StainingSet"), uploadReport)


