library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)
library(limma)

synapseLogin()

reportType <- "Preprocessing"
dataType <- "Quantitative Imaging"
#dataType <- "Metadata"
activityName <- "Annotate data"
#activityName <- "Normalize Cell Data"
#reportDesc <- "Analysis"

## USE TOKEN TO AUTHENTICATE FOR YOUR CURRENT SESSION
setGithubToken("4de3666dc4be8a8174192b1c892202a99311c7ab")
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
                  forceVersion=FALSE)
                  #,
                  #executed=scriptLink)
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
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC240L"), 2),
             StainingSet=c("SS1","SS4"),
             Drug="none",
             Preprocess="av1.6",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC122L"), 2),
             StainingSet=c("SS1","SS4"),
             Drug="none",
             Preprocess="av1.6",
             Segmentation="v2",
             stringsAsFactors=FALSE))

#Get the filepaths for the non-level 1data on the graylab server
filePaths <- grep("Level3",grep("HMEC",dir("../AnnotatedData",full.names = TRUE), value=TRUE),value=TRUE)
#Get annotations from filename
splits <- strsplit2(filePaths,"_")
serverFiles <- data.frame(FilePaths=splits[,1],ss=splits[,2], Level=gsub("Level|.txt","",splits[,6]), stringsAsFactors = FALSE)

  
files <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(files)<-gsub("file.","",colnames(files))

ssDatasets <- merge(ssDatasets,files)
ssDatasets <- select(ssDatasets, -name)
ssDatasets <- reshape2::dcast(ssDatasets,CellLine+StainingSet+Drug+Preprocess+Segmentation~Level, value.var="id")

uploadFile <- function(x){
  dataDir <- paste("../AnnotatedData")
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
                         Level=1,
                         Level1SynID=x[[1]],
                         stringsAsFactors = FALSE)
  
  dataFile$filename <- paste(dataDir,paste0(paste(x$CellLine,x$StainingSet,x$Drug,x$Segmentation,x$Preprocess,"Level1",sep="_"),".txt"),sep="/")
  uploadToSynapse(dataFile, parentId=synapseAnnotatedDataDir)
  
}

res <- dlply(ssDatasets[grepl("HMEC",ssDatasets$CellLine),], c("CellLine","StainingSet"), uploadFile)

#Task
#Upload selected level 3,4 and CS data and reports from graylab server to Synapse
#Get data files available on server
#Use a manual selection list to upload the files
#Create provenance on level 4 and combined files


